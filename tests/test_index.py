import logging

import pytest

import pyfmindex as pfmi

logger = logging.getLogger(__name__)

SEQUENCE = "TACTGTCTTATGAAGATAAGTGAGATAATCTTGACCTGTAGCACTCAGCAGCTGCTGTATTTACCAGGTACAGATAAGACAACA"
MER3 = "CTG"
MER4 = "TACT"
KMERS = [MER3, MER4]
SUFFIX_ARRAY_COMPRESSION_RATIO = 8
KMER_LENGTH_IN_SEED_TABLE = 12
ALPHABET_TYPE = 2


def test_read_index_from_file():
    index = pfmi.read_index_from_file("./tests/index.awfmi", False)
    assert index is not None


def test_find_search_range_for_string(index):
    search_range = index.find_search_range_for_string(MER3)
    assert search_range is not None


def test_create_kmer_search_list():
    kmer_search_list = pfmi.KmerSearchList(5)
    kmer_search_list.fill(KMERS)
    assert kmer_search_list.count == 2


def test_parallel_search_locate(index):
    kmer_search_list = pfmi.KmerSearchList(5)
    kmer_search_list.fill(KMERS)
    kmer_search_list.parallel_search_locate(index)
    for i in range(kmer_search_list.count):
        assert kmer_search_list.kmer_search_data[i].kmer_string.decode() in KMERS


def test_parallel_search_count(index):
    kmer_search_list = pfmi.KmerSearchList(5)
    kmer_search_list.fill(KMERS)
    kmer_search_list.parallel_search_count(index)
    kmers_count = (4, 1)
    for i in range(kmer_search_list.count):
        assert kmer_search_list.kmer_search_data[i].count in kmers_count


def test_read_sequence_from_file(index):
    segment = index.read_sequence_from_file(10, 10)
    assert segment == "TGAAGATAAG"


@pytest.fixture(scope="session")
def config():
    return pfmi.IndexConfiguration(
        SUFFIX_ARRAY_COMPRESSION_RATIO,
        KMER_LENGTH_IN_SEED_TABLE,
        ALPHABET_TYPE,
        True,
        True,
    )


@pytest.fixture(scope="session")
def index(config):
    return pfmi.Index(config, "./tests/index.awfmi", SEQUENCE)
