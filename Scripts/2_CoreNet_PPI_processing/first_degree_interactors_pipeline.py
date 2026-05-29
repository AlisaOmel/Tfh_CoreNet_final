#!/usr/bin/env python
"""
Direct-interactor (DI) permutation pipeline.

The previous version called this "FD" / "first-degree interactors";
the underlying logic is identical, the term is now "DI / direct interactors".

Permutation modes (controlled by --mode):
  di_then_match : find DI of sig_genes once, then build size-matched
                  randoms of the DI set (1000 perms).
  match_then_di : draw a random set the same size as sig_genes, then
                  compute DI of that random set, repeat 1000x.
  sig_only      : skip DI entirely. Compute overlap of sig_genes directly,
                  plus 1000 size-matched random draws.
  both          : run di_then_match + match_then_di (back-compat).
  all           : run all three pipelines.
  none          : skip every permutation pipeline. Useful with --save_di
                  if you only want the real DI set saved.

Independent toggles:
  --save_di         : also save the real DI set of sig_genes (the unique-genes pkl,
                      csv, and protein-pair csv). Runs alongside any --mode value,
                      including --mode none. Note: pipeline_di_then_match no longer
                      saves the DI set itself -- use --save_di if you want it.
  --save_cytoscape  : also save Cytoscape edge + label files for sig_genes' network.
                      Variant chosen by --cytoscape_variant. Independent of --mode.
  --save_formats    : comma-separated formats for gene-list outputs.
                      Default 'csv'. Choices: csv, pkl, txt.

Outputs per pipeline:
  * the per-permutation overlap dataframe
  * the "average overlap" gene list (via get_random_gene_list logic)
  * the real DI set (di_then_match always; --save_di on demand)
  * Cytoscape edge/label/unique-gene files (--save_cytoscape on demand)
"""

import argparse
import csv
import os
import pickle
import random
from concurrent.futures import ProcessPoolExecutor, as_completed

import h5py
import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Original functions (kept as-is, except `random_set` now takes the pool
# explicitly so we don't depend on a module global at call time)
# ---------------------------------------------------------------------------

def make_hint_df(hint_data):
    gene_sets = []
    for entry in hint_data:
        prot_1 = entry[0].decode('utf-8')
        prot_2 = entry[1].decode('utf-8')
        gene_sets.append([prot_1, prot_2])
    return gene_sets


def sig_set_nodot(gene_set, set_name='PPS', output_dir='.', save='True', taiji='False'):
    nodot = []
    if taiji == 'True':
        gene_list = pd.read_csv('{0}'.format(gene_set))
        gene_list = set(gene_list['Genes'])
    else:
        gene_list = list(pd.read_csv('{0}'.format(gene_set)))
        gene_list = set(gene_list)
    for gene in gene_list:
        new_gene = gene.split(".")[0]
        if new_gene not in nodot:
            nodot.append(new_gene)
    if save == 'True':
        with open("{0}/{1}_nodot_genes.pkl".format(output_dir, set_name), "wb") as fp:
            pickle.dump(nodot, fp)
    return nodot


def first_degree_interactors(sig_genes, hint_data, output_dir='.', set_name='PPS',
                             save='off', save_formats=('csv',)):
    """Find direct interactors. When save='on', writes:
      * protein_pairs_df.csv (always, when save='on')
      * unique_genes outputs in formats listed in save_formats
        - 'csv' -> {set_name}_unique_genes_df.csv (single-column 'Genes')
        - 'pkl' -> {set_name}_unique_genes.pkl   (pickled list)
    """
    gene_sets = []
    genes = []
    for entry in hint_data:
        prot_1 = entry[0].decode('utf-8')
        prot_2 = entry[1].decode('utf-8')

        if prot_1 in sig_genes or prot_2 in sig_genes:
            gene_sets.append([prot_1, prot_2])
            if prot_1 not in genes:
                genes.append(prot_1)
            if prot_2 not in genes:
                genes.append(prot_2)

    gene_set_df = pd.DataFrame(gene_sets, columns=['Prot_1', 'Prot_2'])
    my_genes_df = pd.DataFrame(genes, columns=['Genes'])
    if save == 'on':
        print('Saving direct interactors (formats={0})...'.format(tuple(save_formats)))
        # Protein-pair table is always CSV — no other format makes sense for a pair table.
        gene_set_df.to_csv('{0}/{1}_protein_pairs_df.csv'.format(output_dir, set_name),
                           sep=',', index=False)
        if 'csv' in save_formats:
            my_genes_df.to_csv('{0}/{1}_unique_genes_df.csv'.format(output_dir, set_name),
                               sep=',', index=False)
        if 'pkl' in save_formats:
            with open('{0}/{1}_unique_genes.pkl'.format(output_dir, set_name), 'wb') as fp:
                pickle.dump(genes, fp)
    return gene_set_df, genes, my_genes_df


def find_overlap(gene_list, list_name='patternB', output='percent'):
    # gene_list expected as a set
    vinuesa_overlap = vinuesa_list & gene_list
    vinuesa_overlap_percent = len(vinuesa_overlap) / len(vinuesa_list)

    crotty_overlap = crotty_mouse & gene_list
    crotty_overlap_percent = len(crotty_overlap) / len(crotty_mouse)

    hart_overlap = hart_list & gene_list
    hart_overlap_percent = len(hart_overlap) / len(hart_list)

    if output == 'genes':
        return vinuesa_overlap, hart_overlap, crotty_overlap
    if output == 'gene_len':
        return len(vinuesa_overlap), len(hart_overlap), len(crotty_overlap)
    if output == 'genes_gene_len':
        return len(vinuesa_overlap), len(hart_overlap), len(crotty_overlap), gene_list
    if output == 'percent':
        return vinuesa_overlap_percent, hart_overlap_percent, crotty_overlap_percent


def random_set(gene_set, pool):
    """Draw a random set the same size as gene_set from `pool`."""
    random_size_matched_list = set(random.sample(list(pool), len(gene_set)))
    return random_size_matched_list


def _materialize_hint(hint_data):
    """Convert an h5py dataset to a list of (str, str) tuples so it can
    be pickled across process boundaries. Cheap one-time cost; HINT fits
    easily in memory (<1M edges)."""
    out = []
    for entry in hint_data:
        a, b = entry[0], entry[1]
        if isinstance(a, bytes):
            a = a.decode('utf-8')
        if isinstance(b, bytes):
            b = b.decode('utf-8')
        out.append((a, b))
    return out


def _perm_worker(args):
    """Worker for one permutation. Must be top-level so it pickles.
    `hint_pairs` is a list of (str, str) tuples (already materialized from h5py).
    `ref_sets` is a (vinuesa, hart, crotty) tuple of frozensets.
    Returns a 1-row dict in the same shape as the serial loop."""
    (i, sig_genes, pool, do_fd, hint_pairs, ref_sets, set_name, seed) = args

    # Per-task seed -- deterministic given (master_seed, i).
    rng = random.Random(seed)

    # Inline draw to avoid using the global `random` state in workers.
    pool_list = list(pool) if not isinstance(pool, list) else pool
    random_gene_set = set(rng.sample(pool_list, len(sig_genes)))

    if do_fd:
        # Inline DI scan over the materialized list (no h5py needed in the worker).
        sig_set = random_gene_set
        di_genes_list = []
        seen = set()
        for p1, p2 in hint_pairs:
            if p1 in sig_set or p2 in sig_set:
                if p1 not in seen:
                    seen.add(p1); di_genes_list.append(p1)
                if p2 not in seen:
                    seen.add(p2); di_genes_list.append(p2)
        eval_set = seen
    else:
        eval_set = random_gene_set

    vinuesa, hart, crotty = ref_sets
    v = len(vinuesa & eval_set)
    h = len(hart & eval_set)
    c = len(crotty & eval_set)

    return {
        'i': i,
        'Set': 'Random_{0}_{1}'.format(set_name, i),
        'Vinuesa': v,
        'Hart': h,
        'Crotty': c,
        'Genes': eval_set,
    }


def random_overlaps_perm(sig_genes, set_name='PatternB', num_perm=1000,
                         output='genes_gene_len', pool=None, do_fd=False,
                         hint_data=None, n_jobs=1, seed=None):
    """
    Run `num_perm` permutations.
      - do_fd=False: each iteration draws a random set the same size as `sig_genes`
        from `pool` and computes overlaps directly.  (Original behavior.)
      - do_fd=True : each iteration draws a random set the same size as `sig_genes`
        from `pool`, computes its direct interactors (DI), then computes overlaps
        on that DI set.

    Parallelism:
      n_jobs=1  -> serial loop, identical to the original behavior (default).
      n_jobs>1  -> ProcessPoolExecutor across `n_jobs` workers, but ONLY when
                   do_fd=True (the slow regime). For do_fd=False the parallel
                   overhead exceeds the per-task work, so we run serially and
                   print a notice.

    Reproducibility:
      Each task gets seed (master_seed * 1_000_003 + i) so a given (seed, num_perm)
      produces the same draws regardless of n_jobs. Note: results from a parallel
      run will NOT bit-match the original serial run that uses the global RNG --
      use --n_jobs 1 to reproduce previously-published numbers exactly.
    """
    if pool is None:
        raise ValueError("`pool` (background gene universe) must be provided.")
    if do_fd and hint_data is None:
        raise ValueError("`hint_data` must be provided when do_fd=True.")

    use_parallel = n_jobs > 1 and do_fd
    if n_jobs > 1 and not do_fd:
        print('Note: n_jobs={0} requested but do_fd=False; running serially '
              '(parallel overhead > per-task work).'.format(n_jobs))

    if use_parallel:
        # Materialize hint_data once so workers can receive it.
        hint_pairs = _materialize_hint(hint_data)
        ref_sets = (frozenset(vinuesa_list), frozenset(hart_list), frozenset(crotty_mouse))
        master = seed if seed is not None else random.randint(0, 2**31 - 1)
        tasks = [
            (i, set(sig_genes), list(pool), do_fd, hint_pairs, ref_sets,
             set_name, master * 1_000_003 + i)
            for i in range(num_perm)
        ]
        results = []
        with ProcessPoolExecutor(max_workers=n_jobs) as ex:
            for k, res in enumerate(ex.map(_perm_worker, tasks, chunksize=max(1, num_perm // (n_jobs * 4)))):
                results.append(res)
                if k % 50 == 0:
                    print('Num_perm completed (parallel):', k)
        # Sort by i so output order matches the serial version.
        results.sort(key=lambda r: r['i'])
        rows = []
        for r in results:
            rows.append(pd.DataFrame({
                'Set': [r['Set']],
                'Vinuesa': [r['Vinuesa']],
                'Hart': [r['Hart']],
                'Crotty': [r['Crotty']],
                'Genes': [r['Genes']],
            }))
        return pd.concat(rows, axis=0).reset_index(drop=True)

    # ---------- Serial path (unchanged behavior) ----------
    rows = []
    for i in range(0, num_perm):
        print('Num_perm=', i)
        random_gene_set = random_set(sig_genes, pool=pool)

        if do_fd:
            _, di_genes, _ = first_degree_interactors(
                random_gene_set, hint_data,
                set_name='{0}_perm{1}'.format(set_name, i), save='off'
            )
            eval_set = set(di_genes)
        else:
            eval_set = random_gene_set

        random_genes_overlapping = find_overlap(eval_set, output=output)
        df1 = {
            'Set': ['Random_{0}_{1}'.format(set_name, i)],
            'Vinuesa': [random_genes_overlapping[0]],
            'Hart': [random_genes_overlapping[1]],
            'Crotty': [random_genes_overlapping[2]],
            'Genes': [random_genes_overlapping[3]],
        }
        rows.append(pd.DataFrame(df1))

    random_percents_df = pd.concat(rows, axis=0).reset_index(drop=True)
    return random_percents_df


def get_random_gene_list(random_overlap_set, output_dir='.', set_name='PPS',
                         save='on', save_formats=('csv',)):
    random_overlap_set = random_overlap_set.reset_index().drop(columns=['index'])
    random_overlap_subset_human = random_overlap_set.loc[
        (random_overlap_set['Vinuesa'].eq(np.round(np.mean(random_overlap_set['Vinuesa']))))
        & (random_overlap_set['Hart'].eq(np.round(np.mean(random_overlap_set['Hart']))))
    ].reset_index()
    random_overlap_subset_mouse = random_overlap_set.loc[
        (random_overlap_set['Crotty'].eq(np.round(np.mean(random_overlap_set['Crotty']))))
    ].reset_index()

    if len(random_overlap_subset_human) == 0:
        print("WARNING: no permutation matched the rounded human mean overlaps. "
              "Falling back to the row closest to the mean.")
        means = (random_overlap_set[['Vinuesa', 'Hart']].mean())
        dist = ((random_overlap_set[['Vinuesa', 'Hart']] - means) ** 2).sum(axis=1)
        random_overlap_subset_human = random_overlap_set.loc[[dist.idxmin()]].reset_index()
    if len(random_overlap_subset_mouse) == 0:
        print("WARNING: no permutation matched the rounded mouse mean overlaps. "
              "Falling back to the row closest to the mean.")
        means = (random_overlap_set[['Crotty']].mean())
        dist = ((random_overlap_set[['Crotty']] - means) ** 2).sum(axis=1)
        random_overlap_subset_mouse = random_overlap_set.loc[[dist.idxmin()]].reset_index()

    random_overlap_genes_human = random_overlap_subset_human.loc[0]['Genes']
    random_overlap_genes_mouse = random_overlap_subset_mouse.loc[0]['Genes']

    gene_list_human = [g for g in random_overlap_genes_human]
    gene_list_mouse = [g for g in random_overlap_genes_mouse]

    random_overlap_genes_df_human = pd.DataFrame({'Genes': gene_list_human})
    random_overlap_genes_df_mouse = pd.DataFrame({'Genes': gene_list_mouse})

    if save == 'on':
        print('Saving avg-overlap gene lists (formats={0})...'.format(tuple(save_formats)))
        for which, raw, df in (('human', random_overlap_genes_human, random_overlap_genes_df_human),
                               ('mouse', random_overlap_genes_mouse, random_overlap_genes_df_mouse)):
            base = '{0}/{1}_random_upset_{2}'.format(output_dir, set_name, which)
            if 'csv' in save_formats:
                df.to_csv('{0}.csv'.format(base), index=False)
            if 'pkl' in save_formats:
                with open('{0}.pkl'.format(base), 'wb') as fp:
                    pickle.dump(raw, fp)
            if 'txt' in save_formats:
                # one line, comma-separated, no header — legacy format
                with open('{0}.txt'.format(base), 'w', newline='') as f:
                    csv.writer(f).writerow(raw)

    return random_overlap_genes_human, random_overlap_genes_mouse


# ---------------------------------------------------------------------------
# Helpers for input loading
# ---------------------------------------------------------------------------

def _load_gene_csv(path):
    """Load a literature/reference gene set from a CSV file with a single
    column. Tries common column names ('Genes', 'Gene', 'gene', 'symbol'),
    falls back to the first column. Returns a set of strings."""
    df = pd.read_csv(path)
    for c in ('Genes', 'Gene', 'gene', 'genes', 'Symbol', 'symbol'):
        if c in df.columns:
            return set(df[c].dropna().astype(str))
    return set(df.iloc[:, 0].dropna().astype(str))


def load_sig_genes(path, gene_col=None):
    """Load significant gene list from a csv/tsv/pkl/txt file."""
    if path.endswith('.pkl'):
        with open(path, 'rb') as fp:
            obj = pickle.load(fp)
        return set(obj)
    df = pd.read_csv(path)
    if gene_col is not None:
        return set(df[gene_col].dropna().astype(str))
    # try common column names, else use the first column
    for c in ('Genes', 'Gene', 'gene', 'genes', 'Symbol', 'symbol'):
        if c in df.columns:
            return set(df[c].dropna().astype(str))
    return set(df.iloc[:, 0].dropna().astype(str))


# ---------------------------------------------------------------------------
# Cytoscape network builders
# ---------------------------------------------------------------------------

def cytoscape_simple(sig_genes, hint_data):
    """sig_genes + their DI. Returns (edges_df, labels_df, all_genes).
    Labels: 0 = sig_gene, 1 = DI-only (interactor of a sig_gene but not itself sig).
    """
    sig_set = set(sig_genes)
    edges = []
    all_genes = []
    for entry in hint_data:
        prot_1 = entry[0].decode('utf-8')
        prot_2 = entry[1].decode('utf-8')
        if prot_1 in sig_set or prot_2 in sig_set:
            edges.append([prot_1, prot_2])
            if prot_1 not in all_genes:
                all_genes.append(prot_1)
            if prot_2 not in all_genes:
                all_genes.append(prot_2)

    labels = []
    for g in all_genes:
        labels.append([g, 0 if g in sig_set else 1])

    edges_df = pd.DataFrame(edges, columns=['Prot_1', 'Prot_2'])
    labels_df = pd.DataFrame(labels, columns=['Prot', 'Label'])
    return edges_df, labels_df, all_genes


def cytoscape_goldstand(gene_set, sig_list, hint_data):
    """Port of first_degree_goldstand_set. Keep an edge if at least one
    protein is in `gene_set` AND at least one is in `sig_list`.
    Labels: 0 = in both, 1 = sig_list only, 2 = gene_set only, 3 = neither.
    """
    gene_set = set(gene_set)
    sig_list = set(sig_list)
    edges = []
    all_genes = []
    for entry in hint_data:
        prot_1 = entry[0].decode('utf-8')
        prot_2 = entry[1].decode('utf-8')
        if (prot_1 in gene_set or prot_2 in gene_set) \
                and (prot_1 in sig_list or prot_2 in sig_list):
            edges.append([prot_1, prot_2])
            if prot_1 not in all_genes:
                all_genes.append(prot_1)
            if prot_2 not in all_genes:
                all_genes.append(prot_2)

    labels = []
    for g in all_genes:
        if g in gene_set and g in sig_list:
            labels.append([g, 0])
        elif g not in gene_set and g in sig_list:
            labels.append([g, 1])
        elif g in gene_set and g not in sig_list:
            labels.append([g, 2])
        else:
            labels.append([g, 3])

    edges_df = pd.DataFrame(edges, columns=['Prot_1', 'Prot_2'])
    labels_df = pd.DataFrame(labels, columns=['Prot', 'Label'])
    return edges_df, labels_df, all_genes


def cytoscape_set(gene_set, sig_list, sig_di, hint_data):
    """Port of first_degree_cytoscape_set. Keep an edge if both proteins
    are in `gene_set`, OR if exactly one is in `gene_set` AND both are in
    `sig_list`. Labels are computed against `sig_di` (the DI set of sig_genes).
    Labels: 0 = in both gene_set and sig_di, 1 = sig_di only,
            2 = gene_set only, 3 = neither.
    """
    gene_set = set(gene_set)
    sig_list = set(sig_list)
    sig_di = set(sig_di)
    edges = []
    all_genes = []
    for entry in hint_data:
        prot_1 = entry[0].decode('utf-8')
        prot_2 = entry[1].decode('utf-8')

        keep = False
        if prot_1 in gene_set and prot_2 in gene_set:
            keep = True
        elif (prot_1 in gene_set) ^ (prot_2 in gene_set):  # exactly one
            if prot_1 in sig_list and prot_2 in sig_list:
                keep = True

        if keep:
            edges.append([prot_1, prot_2])
            if prot_1 not in all_genes:
                all_genes.append(prot_1)
            if prot_2 not in all_genes:
                all_genes.append(prot_2)

    labels = []
    for g in all_genes:
        if g in gene_set and g in sig_di:
            labels.append([g, 0])
        elif g not in gene_set and g in sig_di:
            labels.append([g, 1])
        elif g in gene_set and g not in sig_di:
            labels.append([g, 2])
        else:
            labels.append([g, 3])

    edges_df = pd.DataFrame(edges, columns=['Prot_1', 'Prot_2'])
    labels_df = pd.DataFrame(labels, columns=['Prot', 'Label'])
    return edges_df, labels_df, all_genes


def _save_cytoscape(edges_df, labels_df, all_genes, prefix, variant, output_dir):
    tag = '{0}_cyto_{1}'.format(prefix, variant)
    edges_df.to_csv('{0}/{1}_edges.csv'.format(output_dir, tag), index=False)
    labels_df.to_csv('{0}/{1}_labels.csv'.format(output_dir, tag), index=False)
    with open('{0}/{1}_unique_genes.pkl'.format(output_dir, tag), 'wb') as fp:
        pickle.dump(all_genes, fp)
    print('  Cytoscape {0}: {1} edges, {2} unique genes'.format(
        variant, len(edges_df), len(all_genes)))


# ---------------------------------------------------------------------------
# Pipelines
# ---------------------------------------------------------------------------

def pipeline_save_di(sig_genes, hint_data, prefix, output_dir, save_formats=('csv',)):
    """Just compute and save the real DI set of sig_genes.
    No permutations, no overlap calculation. Useful with --mode none
    when you only want the DI set on disk."""
    print('\n=== Saving real DI (direct interactor) set of sig_genes ===')
    set_name = '{0}_real_DI'.format(prefix)

    di_pairs_df, di_genes, di_genes_df = first_degree_interactors(
        sig_genes, hint_data, output_dir=output_dir,
        set_name=set_name, save='on', save_formats=save_formats
    )
    print('Real DI set size: {0}'.format(len(di_genes)))
    return di_genes


def pipeline_di_then_match(sig_genes, hint_data, all_hint_genes, prefix,
                           output_dir, num_perm, save_formats=('csv',),
                           n_jobs=1, seed=None):
    """First find DI of sig_genes, then size-match-randomize.
    NOTE: the DI set itself is NOT saved here -- use TRUE_DI / pipeline_save_di
    if you want it on disk. We only compute it in memory for sizing the perm."""
    print('\n=== Pipeline 1: first_DI_then_size_match ===')
    set_name = '{0}_first_DI_then_size_match'.format(prefix)

    # 1. DI of sig_genes (in-memory only)
    di_pairs_df, di_genes, di_genes_df = first_degree_interactors(
        sig_genes, hint_data, output_dir=output_dir,
        set_name=set_name, save='off'
    )
    print('DI set size: {0}'.format(len(di_genes)))

    # 2. random size-matched sets matching the DI set size
    # do_fd=False here, so n_jobs has no effect (parallel overhead > work).
    perm_df = random_overlaps_perm(
        di_genes, set_name=set_name, num_perm=num_perm,
        output='genes_gene_len', pool=all_hint_genes, do_fd=False,
        n_jobs=n_jobs, seed=seed,
    )
    perm_df.to_csv('{0}/{1}_random_overlaps_df.csv'.format(output_dir, set_name), index=False)

    # 3. avg-overlap gene list
    get_random_gene_list(perm_df, output_dir=output_dir, set_name=set_name,
                         save='on', save_formats=save_formats)


def pipeline_match_then_di(sig_genes, hint_data, all_hint_genes, prefix,
                           output_dir, num_perm, save_formats=('csv',),
                           n_jobs=1, seed=None):
    """Draw size-matched random of sig_genes first, then find DI per perm.
    This is the slow pipeline -- n_jobs > 1 actually helps here."""
    print('\n=== Pipeline 2: first_size_match_then_DI ===')
    set_name = '{0}_first_size_match_then_DI'.format(prefix)

    perm_df = random_overlaps_perm(
        sig_genes, set_name=set_name, num_perm=num_perm,
        output='genes_gene_len', pool=all_hint_genes,
        do_fd=True, hint_data=hint_data,
        n_jobs=n_jobs, seed=seed,
    )
    perm_df.to_csv('{0}/{1}_random_overlaps_df.csv'.format(output_dir, set_name), index=False)

    get_random_gene_list(perm_df, output_dir=output_dir, set_name=set_name,
                         save='on', save_formats=save_formats)


def pipeline_sig_only(sig_genes, all_hint_genes, prefix, output_dir, num_perm,
                      save_formats=('csv',), n_jobs=1, seed=None):
    """No DI step on either side. Draw randoms at len(sig_genes) from the
    background pool. Real sig_genes overlap is computed elsewhere."""
    print('\n=== Pipeline 3: sig_only_size_match ===')
    set_name = '{0}_sig_only_size_match'.format(prefix)

    perm_df = random_overlaps_perm(
        sig_genes, set_name=set_name, num_perm=num_perm,
        output='genes_gene_len', pool=all_hint_genes, do_fd=False,
        n_jobs=n_jobs, seed=seed,
    )
    perm_df.to_csv('{0}/{1}_random_overlaps_df.csv'.format(output_dir, set_name), index=False)

    get_random_gene_list(perm_df, output_dir=output_dir, set_name=set_name,
                         save='on', save_formats=save_formats)


def pipeline_save_cytoscape(sig_genes, hint_data, prefix, output_dir,
                            variant='simple', gold_standard=None):
    """Build Cytoscape-format edge + label files.

    variant:
      'simple'    -- sig_genes + their DI. 2 labels.
      'goldstand' -- requires gold_standard. Edges where any node is in
                     gold_standard AND any node is in sig_genes. 4 labels.
      'cytoscape' -- requires gold_standard. Edges within gold_standard
                     OR exactly-one-in edges where both nodes are in sig_genes.
                     Labels computed against the DI of sig_genes. 4 labels.
      'all'       -- run all three variants.
    """
    print('\n=== Building Cytoscape network(s): variant={0} ==='.format(variant))

    variants_to_run = ['simple', 'goldstand', 'cytoscape'] if variant == 'all' else [variant]

    needs_gold = any(v in ('goldstand', 'cytoscape') for v in variants_to_run)
    if needs_gold and gold_standard is None:
        raise ValueError(
            "variants 'goldstand' and 'cytoscape' require `gold_standard`.")

    # The 'cytoscape' variant needs the DI of sig_genes for labeling.
    # Compute fresh here so we don't depend on call order with other pipelines.
    sig_di = None
    if 'cytoscape' in variants_to_run:
        _, sig_di, _ = first_degree_interactors(sig_genes, hint_data, save='off')

    for v in variants_to_run:
        if v == 'simple':
            edges_df, labels_df, genes = cytoscape_simple(sig_genes, hint_data)
        elif v == 'goldstand':
            edges_df, labels_df, genes = cytoscape_goldstand(
                gold_standard, sig_genes, hint_data
            )
        elif v == 'cytoscape':
            edges_df, labels_df, genes = cytoscape_set(
                gold_standard, sig_genes, sig_di, hint_data
            )
        else:
            raise ValueError("unknown cytoscape variant: {0}".format(v))
        _save_cytoscape(edges_df, labels_df, genes, prefix, v, output_dir)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--sig_genes', required=True,
                   help='Path to significant genes file (.csv/.tsv/.txt/.pkl).')
    p.add_argument('--prefix', required=True,
                   help='Prefix used in all output filenames (e.g. PPS, taiji).')
    p.add_argument('--gene_col', default=None,
                   help='Column name holding gene symbols. If omitted, common names are tried.')
    p.add_argument('--hint_h5', required=True,
                   help='Path to HINT interactions .h5 file.')
    p.add_argument('--all_hint_genes', required=True,
                   help='Path to pickled list of all HINT background genes.')
    p.add_argument('--vinuesa', required=True,
                   help='CSV of Vinuesa gene set (single column, header e.g. "Genes").')
    p.add_argument('--hart', required=True,
                   help='CSV of Hart gene set (single column, header e.g. "Genes").')
    p.add_argument('--crotty', required=True,
                   help='CSV of Crotty (mouse) gene set (single column, header e.g. "Genes").')
    p.add_argument('--output_dir', required=True, help='Where to save outputs.')
    p.add_argument('--num_perm', type=int, default=1000, help='Number of permutations.')
    p.add_argument('--mode',
                   choices=['di_then_match', 'match_then_di', 'sig_only',
                            'both', 'all', 'none'],
                   default='both',
                   help="Which permutation pipeline(s) to run. "
                        "'both' = the two DI pipelines (back-compat). "
                        "'all' = all three pipelines including sig_only. "
                        "'none' = skip every permutation pipeline (use with "
                        "--save_di to only save the real DI set).")
    p.add_argument('--save_di', action='store_true',
                   help='Independent toggle: also save the real DI (direct '
                        'interactor) set of sig_genes. Runs alongside any '
                        '--mode value, including --mode none.')
    p.add_argument('--save_cytoscape', action='store_true',
                   help='Independent toggle: also save Cytoscape edge + label '
                        'files for the sig_genes network. Variant set by '
                        '--cytoscape_variant.')
    p.add_argument('--cytoscape_variant',
                   choices=['simple', 'goldstand', 'cytoscape', 'all'],
                   default='simple',
                   help="Which Cytoscape network to build. "
                        "'simple' = sig + DI (no gold standard needed). "
                        "'goldstand' / 'cytoscape' = use gold-standard list. "
                        "'all' = build all three.")
    p.add_argument('--gold_standard', default=None,
                   help='CSV of gold-standard gene set (single column, header '
                        'e.g. "Genes"). Required for cytoscape_variant in '
                        "{'goldstand', 'cytoscape', 'all'}.")
    p.add_argument('--seed', type=int, default=42,
                   help='Master random seed. Each iteration advances the RNG, '
                        'so 1000 perms produce 1000 different draws.')
    p.add_argument('--save_formats', default='csv',
                   help="Comma-separated list of formats for saving gene lists. "
                        "Choices: csv, pkl, txt. Default: 'csv' (only the format "
                        "the downstream R scripts read). Example: --save_formats csv,pkl")
    p.add_argument('--n_jobs', type=int, default=1,
                   help="Parallel workers for permutations. Default 1 (serial). "
                        "Only speeds up pipeline 2 (match_then_di) where each perm "
                        "runs a full HINT scan. Pipelines 1 and 3 ignore this.")
    args = p.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    random.seed(args.seed)

    # --- load reference sets as module-level globals so find_overlap sees them
    # Literature sets are now CSV files with a single 'Genes' column.
    global vinuesa_list, hart_list, crotty_mouse
    vinuesa_list = _load_gene_csv(args.vinuesa)
    hart_list = _load_gene_csv(args.hart)
    crotty_mouse = _load_gene_csv(args.crotty)

    # --- HINT
    ppi = h5py.File(args.hint_h5, 'r')
    hint_data = ppi['edges']
    with open(args.all_hint_genes, 'rb') as fp:
        all_hint_genes = pickle.load(fp)

    # --- input
    sig_genes = load_sig_genes(args.sig_genes, gene_col=args.gene_col)
    print('Loaded {0} significant genes.'.format(len(sig_genes)))

    # --- run
    save_formats = tuple(s.strip() for s in args.save_formats.split(',') if s.strip())
    valid = {'csv', 'pkl', 'txt'}
    bad = [f for f in save_formats if f not in valid]
    if bad:
        raise ValueError("Unknown save_formats {0}; valid: {1}".format(bad, sorted(valid)))

    if args.save_di:
        pipeline_save_di(sig_genes, hint_data, args.prefix, args.output_dir,
                         save_formats=save_formats)

    if args.save_cytoscape:
        gold = None
        if args.gold_standard is not None:
            gold = _load_gene_csv(args.gold_standard)
        pipeline_save_cytoscape(sig_genes, hint_data, args.prefix,
                                args.output_dir,
                                variant=args.cytoscape_variant,
                                gold_standard=gold)

    if args.mode in ('di_then_match', 'both', 'all'):
        pipeline_di_then_match(sig_genes, hint_data, all_hint_genes,
                               args.prefix, args.output_dir, args.num_perm,
                               save_formats=save_formats,
                               n_jobs=args.n_jobs, seed=args.seed)
    if args.mode in ('match_then_di', 'both', 'all'):
        pipeline_match_then_di(sig_genes, hint_data, all_hint_genes,
                               args.prefix, args.output_dir, args.num_perm,
                               save_formats=save_formats,
                               n_jobs=args.n_jobs, seed=args.seed)
    if args.mode in ('sig_only', 'all'):
        pipeline_sig_only(sig_genes, all_hint_genes,
                          args.prefix, args.output_dir, args.num_perm,
                          save_formats=save_formats,
                          n_jobs=args.n_jobs, seed=args.seed)

    if args.mode == 'none' and not args.save_di and not args.save_cytoscape:
        print('\nNothing to do: --mode is "none" and no --save_* flag was set.')

    print('\nDone.')


if __name__ == '__main__':
    main()
