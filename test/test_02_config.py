import os
from gimmemotifs.config import get_build_dir, MotifConfig, parse_denovo_params


def test_get_build_dir():
    build_dir = get_build_dir()
    assert build_dir is None or "gimmemotifs" in build_dir


def test_create_default_config():
    config = MotifConfig()
    os.remove(config.user_config)

    config.create_default_config()
    mbin = config.config.get("MEME", "bin")
    assert os.path.isfile(mbin) and os.access(mbin, os.X_OK)
    assert mbin.endswith("meme")


def test__upgrade_config():
    config = MotifConfig()

    # outdated param. should not exist
    width = config.config.get("params", "width", fallback=None)
    assert width is None

    # add the outdated param
    size = config.config["params"]["size"]
    config.config.set(
        option="width",
        section="params",
        value=size,
    )
    del config.config["params"]["size"]
    width = config.config.get("params", "width", fallback=None)
    assert width == size

    # after the update, width was removed and its value stored as "size"
    config._upgrade_config()
    width = config.config.get("params", "width", fallback=None)
    assert width is None
    newsize = config.config.get("params", "size", fallback=None)
    assert newsize == size


def test_bin_dir():
    config = MotifConfig()

    mbin = config.bin("MEME")
    assert os.path.isfile(mbin) and os.access(mbin, os.X_OK)
    assert mbin.endswith("meme")

    mdir = config.dir("MEME")
    assert os.path.isdir(mdir)
    assert mdir in mbin


def test_list_installed_libraries():
    config = MotifConfig()
    libs = config.list_installed_libraries()
    pfms = [f for f in os.listdir("data/motif_databases") if f.endswith(".pfm")]
    assert libs == sorted(pfms)


def test_parse_denovo_params():
    params = parse_denovo_params()
    assert isinstance(params["background"], list)
    assert params["max_time"] == -1
