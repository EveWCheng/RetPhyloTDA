import numpy as np
import pytest

import sim_main


def _simulate(seed):
    """Run a few small seeded simulations mixing all three hybridization types."""
    np.random.seed(seed)
    return sim_main.sim_bdh_age(
        age=10,
        numbsim=10,
        lambda_=0.1,
        mu=0.1,
        nu=0.3,
        hybprops=[1, 1, 1],
        hyb_inher_fxn=sim_main.hyb_inher_fxn,
        stopping_num_leaves=6,
    )


# gudhi only reads the lower triangle, so an asymmetric matrix would silently corrupt the filtration
@pytest.mark.parametrize("midpoint_top_percent", [None, 10])
def test_build_dist_matrix_is_symmetric_for_simulated_networks(tmp_path, midpoint_top_percent):
    results = _simulate(seed=42)
    assert len(results) > 0

    for i, r in enumerate(results):
        output_path = tmp_path / f"sim{i}"
        output_path.mkdir()

        dist_matrix, _ = sim_main.build_dist_matrix(
            r["tree_only_G"].to_undirected(),
            r["distance_include_tips"],
            r["distance_only_tips"],
            use_data_prep=False,
            weight_attr="length",
            output_path=str(output_path),
            distance_include_tips_changes=r["distance_include_tips_changes"],
            midpoint_top_percent=midpoint_top_percent,
        )

        assert np.allclose(dist_matrix, dist_matrix.T), f"sim{i} distance matrix is not symmetric"


# distance_include_tips should carry exactly the tip-to-tip values from distance_only_tips
def test_distance_only_tips_agrees_with_distance_include_tips():
    results = _simulate(seed=42)
    assert len(results) > 0

    for i, r in enumerate(results):
        for a, row in r["distance_only_tips"].items():
            for b, d in row.items():
                assert r["distance_include_tips"][a][b] == pytest.approx(d), f"sim{i} tips {a}, {b} disagree"
