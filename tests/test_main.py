import cProfile
import os
import pstats
from pstats import SortKey

import main


# Not a correctness check -- run explicitly (e.g. `pytest tests/test_main.py -k profile`)
# to profile a full main() run and inspect cumulative time per function.
def test_profile_main():
    pr = cProfile.Profile()
    pr.enable()
    main.main()
    pr.disable()

    profile_path = os.path.join(main.HERE, "profile")
    with open(profile_path, "w") as f:
        pstats.Stats(pr, stream=f).sort_stats(SortKey.CUMULATIVE).print_stats()
