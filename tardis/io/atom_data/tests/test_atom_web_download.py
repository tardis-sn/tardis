from tardis.io.atom_data.atom_web_download import get_atomic_repo_config


def test_default_atomic_data_uses_raw_github_url():
    atomic_data_name = "kurucz_cd23_chianti_H_He_latest"
    expected_url = (
        "https://raw.githubusercontent.com/tardis-sn/"
        "tardis-regression-data/main/atom_data/"
        f"{atomic_data_name}.h5"
    )

    atomic_repo = get_atomic_repo_config()

    assert atomic_repo[atomic_data_name]["url"] == expected_url
    assert atomic_repo[atomic_data_name]["mirrors"] == [expected_url]
