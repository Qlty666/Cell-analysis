"""Bounded task submission for large ligand libraries."""
from concurrent.futures import FIRST_COMPLETED, wait


def bounded_futures(pool, function, items, max_pending):
    if max_pending < 1:
        raise ValueError("max_pending must be positive")
    iterator = iter(items)
    pending = {}
    exhausted = False
    while pending or not exhausted:
        while not exhausted and len(pending) < max_pending:
            try:
                item = next(iterator)
            except StopIteration:
                exhausted = True
            else:
                pending[pool.submit(function, item)] = item
        if pending:
            done, _ = wait(pending, return_when=FIRST_COMPLETED)
            for future in done:
                yield future, pending.pop(future)
