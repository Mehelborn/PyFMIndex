import logging

import pytest

import pyfmindex as pfmi
from pyfmindex import _pyfmindex as _pfmi

logger = logging.getLogger(__name__)

SEQUENCE = "TACTGTCTTATGAAGATAAGTGAGATAATCTTGACCTGTAGCACTCAGCAGCTGCTGTATTTACCAGGTACAGATAAGACAACA"
MER3 = "CTG"
MER4 = "TACT"
SUFFIX_ARRAY_COMPRESSION_RATIO = 8
KMER_LENGTH_IN_SEED_TABLE = 12
ALPHABET_TYPE = 2


def test_index_config_creation(config):
    assert config is not None


def test_index_creation(config):
    index = pfmi.Index(config, "./tests/index.awfmi", SEQUENCE)
    assert index is not None


def test_read_index_from_file():
    index = pfmi.read_index_from_file("./tests/index.awfmi", False)
    assert index is not None


def test_convert_search_range():
    search_range = _pfmi._SearchRange(5, 10)
    fields = {}
    for field, _ in search_range._fields_:
        fields[field] = getattr(search_range, field)
    search_range = pfmi.SearchRange(**fields)
    assert search_range is not None


def test_find_search_range_for_string(config):
    index = pfmi.Index(config, "./tests/index.awfmi", SEQUENCE)
    search_range = index.find_search_range_for_string(MER3)
    logger.info(search_range)
    assert search_range is not None


def test_create_kmer_search_list(config):
    kmer_search_list: pfmi.KmerSearchList = pfmi.KmerSearchList(5)
    kmer_search_list.fill_out_list(["CTG", "AAT"], [3, 3])
    assert kmer_search_list.count == 2


def test_parallel_search_locate(config):
    index = pfmi.Index(config, "./tests/index.awfmi", SEQUENCE)
    kmer_search_list: pfmi.KmerSearchList = pfmi.KmerSearchList(5)
    kmer_search_list.fill_out_list(["CTG", "AAT"], [3, 3])
    pfmi.parallel_search_locate(index, kmer_search_list, 2)
    ksd = kmer_search_list._kmer_search_list.contents.kmer_search_data[0]
    for i in range(ksd.count):
        logger.info(ksd.position_list[i])
        logger.info(ksd.kmer_string)


def test_read_sequence_from_file(config):
    index = pfmi.Index(config, "./tests/index.awfmi", SEQUENCE)
    segment = index.read_sequence_from_file(10, 10)
    assert segment == "TGAAGATAAG"


@pytest.fixture
def config():
    return pfmi.IndexConfiguration(
        SUFFIX_ARRAY_COMPRESSION_RATIO,
        KMER_LENGTH_IN_SEED_TABLE,
        ALPHABET_TYPE,
        True,
        True,
    )
