import numpy as np

# Debug version of cumulative integrate by blocks
def cumulative_integrate_array_by_blocks_debug(f, x, block_references, normalize=True):
    integrated = np.zeros_like(f)
    dx = np.diff(x)

    for i in range(f.shape[1]):
        # cumulative contributions
        contribs = dx * 0.5 * (f[1:, i] + f[:-1, i])
        cumsum = np.zeros(f.shape[0])
        cumsum[1:] = np.cumsum(contribs)

        for j in range(len(block_references) - 1):
            start = block_references[j]
            stop = block_references[j + 1]
            block = cumsum[start:stop] - cumsum[start]

            # Handle normalization
            if normalize and len(block) > 1 and block[-1] != 0:
                block = block / block[-1]
            elif len(block) == 1:
                block = np.array([f[start, i]])  # preserve original value

            # Debug info
            print(
                f"[DEBUG] start={start}, stop={stop}, "
                f"block.shape={block.shape}, "
                f"target.shape={integrated[start:stop, i].shape}"
            )

            # Assign block into integrated array
            integrated[start:stop, i] = block

    return integrated


if __name__ == "__main__":
    # Dummy input
    f = np.array([
        [1.0, 2.0],
        [2.0, 3.0],
        [3.0, 4.0],
        [4.0, 5.0],
        [5.0, 6.0],
    ])
    x = np.linspace(0, 1, f.shape[0])
    block_refs = np.array([0, 2, 4])  # last block now includes last element

    print("Running debug cumulative integration...")
    result = cumulative_integrate_array_by_blocks_debug(f, x, block_refs, normalize=True)
    print("Result:\n", result)
