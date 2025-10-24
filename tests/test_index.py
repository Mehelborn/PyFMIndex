import logging

import pytest

from pyfmindex import IndexConfiguration, Index, read_index_from_file

logger = logging.getLogger(__name__)

SEQUENCE = "ACGTACGTTTACAGT"
SUFFIX_ARRAY_COMPRESSION_RATIO = 8
KMER_LENGTHIN_SEED_TABLE = 12
ALPHABET_TYPE = 2


def test_index_config_creation(config):
    assert config is not None


def test_index_creation(config):
    index = Index(config, "./tests/index.awfmi", SEQUENCE)
    assert index is not None


def test_read_index_from_file():
    index = read_index_from_file("./tests/index.awfmi", False)
    assert index is not None


@pytest.fixture
def config():
    return IndexConfiguration(
        SUFFIX_ARRAY_COMPRESSION_RATIO,
        KMER_LENGTHIN_SEED_TABLE,
        ALPHABET_TYPE,
        True,
        True,
    )
