"""
Unit tests for PaperCollectionService.
"""

import io
import pytest
from unittest.mock import Mock, patch
from datetime import datetime
from xml.etree import ElementTree as ET

from src.pubmed_miner.services.paper_collection import PaperCollectionService
from src.pubmed_miner.utils.error_handler import APIError


def _parse_esearch_xml(xml_bytes):
    """Parse Entrez esearch XML bytes into dict."""
    root = ET.fromstring(xml_bytes)
    id_list = root.find("IdList")
    if id_list is None:
        return {"IdList": []}
    ids = [eid.text for eid in id_list.findall("Id")]
    return {"IdList": ids}


def _parse_efetch_xml(xml_bytes):
    """Parse Entrez efetch XML bytes returned via Entrez.read() into Paper-like dicts.

    Entrez.read() on efetch XML returns {'PubmedArticle': [dict, ...]} where each
    dict has 'MedlineCitation' -> {'PMID', 'Article': {...}}. This mirrors the
    structure the service's _parse_paper_record expects.
    """
    root = ET.fromstring(xml_bytes)
    articles = []
    for article_elem in root.findall(".//PubmedArticle"):
        medline = article_elem.find("MedlineCitation")
        if medline is None:
            continue
        pmid = medline.findtext("PMID", default="")
        article = medline.find("Article")
        title = article.findtext("ArticleTitle", default="") if article is not None else ""
        authors = []
        author_list = article.find("AuthorList") if article is not None else None
        if author_list is not None:
            for author in author_list.findall("Author"):
                last = author.findtext("LastName", default="")
                fore = author.findtext("ForeName", default="")
                if last and fore:
                    authors.append(f"{fore} {last}")
                elif last:
                    authors.append(last)
        journal = article.findtext("Journal/Title", default="") if article is not None else ""
        date_elem = article.find("ArticleDate") if article is not None else None
        year = int(date_elem.findtext("Year", "2000")) if date_elem is not None else 2000
        month = int(date_elem.findtext("Month", "1")) if date_elem is not None else 1
        day = int(date_elem.findtext("Day", "1")) if date_elem is not None else 1
        from datetime import datetime as dt
        pub_date = dt(year, month, day)
        abstract = ""
        abstract_elem = article.find("Abstract/AbstractText") if article is not None else None
        if abstract_elem is not None:
            abstract = abstract_elem.text or ""
        doi = ""
        pubmed_data = article_elem.find("PubmedData")
        if pubmed_data is not None:
            for aid in pubmed_data.findall(".//ArticleId"):
                if aid.get("IdType") == "doi":
                    doi = aid.text or ""
                    break
        # Build nested dict matching Entrez.read() output structure
        articles.append({
            "MedlineCitation": {
                "PMID": pmid,
                "Article": {
                    "ArticleTitle": title,
                    "AuthorList": [{"LastName": a.split()[-1] if " " in a else a,
                                     "ForeName": " ".join(a.split()[:-1]) if " " in a else ""}
                                    if a else {} for a in authors] if authors else [],
                    "Journal": {"Title": journal},
                    "ArticleDate": {"Year": str(year), "Month": str(month), "Day": str(day)},
                    "Abstract": {"AbstractText": abstract},
                },
            },
            "PubmedData": {
                "ArticleIdList": [
                    {"IdType": "doi", "value": doi}
                ] if doi else []
            },
        })
    return {"PubmedArticle": articles}


def _make_esearch_handle(pmids):
    """Create mock Entrez esearch handle."""
    xml = (
        '<?xml version="1.0" ?>\n'
        '<!DOCTYPE eSearchResult PUBLIC "-//NLM//NLM E-utilities customize//EN" '
        '"https://www.ncbi.nlm.nih.gov/entrez/query/static/eSearchResult.dtd">\n'
        "<eSearchResult>"
        "<IdList>"
        + "".join(f"<Id>{pmid}</Id>" for pmid in pmids)
        + "</IdList>"
        "</eSearchResult>"
    )
    return Mock(
        read=Mock(return_value=xml.encode()),
        close=Mock(),
        __enter__=Mock(return_value=Mock()),
        __exit__=Mock(return_value=False),
    )


def _make_efetch_handle(xml_body):
    """Create mock Entrez efetch handle."""
    xml = (
        '<?xml version="1.0" ?>\n'
        '<!DOCTYPE PubmedArticleSet PUBLIC "-//NLM//NLM PubMed/DTD PubMed Markup 1.0//EN" '
        '"https://www.ncbi.nlm.nih.gov/entrez/query/static/PubmedArticle.dtd">\n'
        + xml_body
    )
    return Mock(
        read=Mock(return_value=xml.encode()),
        close=Mock(),
        __enter__=Mock(return_value=Mock()),
        __exit__=Mock(return_value=False),
    )


class TestPaperCollectionService:
    """Test cases for PaperCollectionService class."""

    def setup_method(self):
        """Set up test fixtures."""
        self.service = PaperCollectionService(
            email="test@example.com",
            rate_limit=10.0,
        )

    def test_initialization(self):
        """Test service initialization."""
        assert self.service.email == "test@example.com"
        assert self.service.rate_limit == 10.0

    @patch("src.pubmed_miner.services.paper_collection.Entrez.read")
    @patch("src.pubmed_miner.services.paper_collection.Entrez.esearch")
    def test_search_papers_success(self, mock_esearch, mock_read):
        """Test successful paper search."""

        def side_effect(handle):
            return _parse_esearch_xml(handle.read())

        mock_read.side_effect = side_effect
        mock_esearch.return_value = _make_esearch_handle(["12345", "67890", "11111"])

        pmids = self.service.search_papers("machine learning", max_results=100)

        assert len(pmids) == 3
        assert "12345" in pmids
        assert "67890" in pmids
        assert "11111" in pmids

        mock_esearch.assert_called_once()
        call_args = mock_esearch.call_args[1]
        assert call_args["db"] == "pubmed"

    @patch("src.pubmed_miner.services.paper_collection.Entrez.esearch")
    def test_search_papers_empty_query(self, mock_esearch):
        """Test search with empty query."""
        with pytest.raises(ValueError, match="Query cannot be empty"):
            self.service.search_papers("", max_results=100)

        with pytest.raises(ValueError, match="Query cannot be empty"):
            self.service.search_papers("   ", max_results=100)

    def test_search_papers_invalid_max_results(self):
        """Test search with invalid max_results."""
        with pytest.raises(ValueError, match="Max results must be positive"):
            self.service.search_papers("test", max_results=0)

        with pytest.raises(ValueError, match="Max results must be positive"):
            self.service.search_papers("test", max_results=-10)

    @patch("src.pubmed_miner.services.paper_collection.Entrez.read")
    @patch("src.pubmed_miner.services.paper_collection.Entrez.esearch")
    def test_search_papers_no_results(self, mock_esearch, mock_read):
        """Test search with no results."""

        def side_effect(handle):
            return _parse_esearch_xml(handle.read())

        mock_read.side_effect = side_effect
        mock_esearch.return_value = _make_esearch_handle([])

        pmids = self.service.search_papers("nonexistent query")

        assert pmids == []

    @patch("src.pubmed_miner.services.paper_collection.Entrez.esearch")
    def test_search_papers_api_error(self, mock_esearch):
        """Test search with API error."""
        mock_esearch.side_effect = Exception("PubMed API error")

        with pytest.raises(APIError, match="PubMed search failed"):
            self.service.search_papers("test query")

    @patch("src.pubmed_miner.services.paper_collection.Entrez.read")
    @patch("src.pubmed_miner.services.paper_collection.Entrez.efetch")
    def test_get_paper_details_success(self, mock_efetch, mock_read):
        """Test successful paper details retrieval."""
        xml_body = (
            "<PubmedArticleSet>"
            "<PubmedArticle>"
            "<MedlineCitation>"
            "<PMID>12345</PMID>"
            "<Article>"
            "<ArticleTitle>Test Paper Title</ArticleTitle>"
            "<AuthorList>"
            "<Author>"
            "<LastName>Doe</LastName>"
            "<ForeName>John</ForeName>"
            "</Author>"
            "<Author>"
            "<LastName>Smith</LastName>"
            "<ForeName>Jane</ForeName>"
            "</Author>"
            "</AuthorList>"
            "<Journal>"
            "<Title>Nature Medicine</Title>"
            "</Journal>"
            "<ArticleDate>"
            "<Year>2023</Year>"
            "<Month>01</Month>"
            "<Day>15</Day>"
            "</ArticleDate>"
            "<Abstract>"
            "<AbstractText>This is a test abstract.</AbstractText>"
            "</Abstract>"
            "</Article>"
            "</MedlineCitation>"
            "<PubmedData>"
            "<ArticleIdList>"
            "<ArticleId IdType=\"doi\">10.1038/test.2023.12345</ArticleId>"
            "</ArticleIdList>"
            "</PubmedData>"
            "</PubmedArticle>"
            "</PubmedArticleSet>"
        )

        def side_effect(handle):
            return _parse_efetch_xml(handle.read())

        mock_read.side_effect = side_effect
        mock_efetch.return_value = _make_efetch_handle(xml_body)

        papers = self.service.get_paper_details(["12345"])

        assert len(papers) == 1
        paper = papers[0]
        assert paper.pmid == "12345"
        assert paper.title == "Test Paper Title"
        assert len(paper.authors) == 2
        assert "John Doe" in paper.authors
        assert "Jane Smith" in paper.authors
        assert paper.journal == "Nature Medicine"
        assert paper.abstract == "This is a test abstract."
        assert paper.doi == "10.1038/test.2023.12345"

    def test_get_paper_details_empty_pmids(self):
        """Test paper details with empty PMID list."""
        papers = self.service.get_paper_details([])
        assert papers == []

    def test_get_paper_details_invalid_pmids(self):
        """Test paper details with invalid PMIDs."""
        with pytest.raises(ValueError, match="All PMIDs must be non-empty strings"):
            self.service.get_paper_details(["12345", "", "67890"])

        with pytest.raises(ValueError, match="All PMIDs must be non-empty strings"):
            self.service.get_paper_details([None, "12345"])

    @patch("src.pubmed_miner.services.paper_collection.Entrez.efetch")
    def test_get_paper_details_api_error(self, mock_efetch):
        """Test paper details with API error."""
        mock_efetch.side_effect = Exception("PubMed API error")

        with pytest.raises(APIError, match="Failed to fetch paper details"):
            self.service.get_paper_details(["12345"])

    @patch("src.pubmed_miner.services.paper_collection.time.sleep")
    def test_rate_limiting(self, mock_sleep):
        """Test rate limiting functionality."""
        service = PaperCollectionService(rate_limit=1.0)

        with patch.object(service, "search_papers") as mock_request:
            mock_request.return_value = []

            service.search_papers("test1")
            service.search_papers("test2")

            assert mock_sleep.call_count >= 1

    def test_parse_authors_various_formats(self):
        """Test author parsing with various formats."""
        authors_xml = """
        <AuthorList>
            <Author>
                <LastName>Doe</LastName>
                <ForeName>John A</ForeName>
            </Author>
            <Author>
                <LastName>Smith</LastName>
                <ForeName>Jane</ForeName>
            </Author>
        </AuthorList>"""

        authors = self.service._parse_authors(authors_xml)
        assert "John A Doe" in authors
        assert "Jane Smith" in authors

    def test_parse_date_various_formats(self):
        """Test date parsing with various formats."""
        date_xml = """
        <ArticleDate>
            <Year>2023</Year>
            <Month>01</Month>
            <Day>15</Day>
        </ArticleDate>"""

        date = self.service._parse_date(date_xml)
        assert date == datetime(2023, 1, 15)

        date_xml_year = """
        <PubDate>
            <Year>2023</Year>
        </PubDate>"""

        date = self.service._parse_date(date_xml_year)
        assert date.year == 2023
        assert date.month == 1
        assert date.day == 1

    def test_batch_processing(self):
        """Test batch processing of PMIDs."""
        pmids = [str(i) for i in range(250)]

        with patch.object(self.service, "_fetch_paper_batch") as mock_fetch:
            mock_fetch.return_value = []

            self.service.get_paper_details(pmids)

            assert mock_fetch.call_count > 1

    @patch("src.pubmed_miner.services.paper_collection.Entrez.read")
    @patch("src.pubmed_miner.services.paper_collection.Entrez.esearch")
    def test_retry_mechanism(self, mock_esearch, mock_read):
        """Test retry mechanism for failed requests."""

        def side_effect(handle):
            return _parse_esearch_xml(handle.read())

        mock_read.side_effect = side_effect
        mock_esearch.side_effect = [
            Exception("Temporary error"),
            _make_esearch_handle([]),
        ]

        pmids = self.service.search_papers("test query")
        assert pmids == []
        assert mock_esearch.call_count == 2

    def test_validate_pmid_format(self):
        """Test PMID format validation."""
        assert self.service._validate_pmid("12345") is True
        assert self.service._validate_pmid("1") is True
        assert self.service._validate_pmid("12345678") is True

        assert self.service._validate_pmid("") is False
        assert self.service._validate_pmid("abc123") is False
        assert self.service._validate_pmid("123456789") is False
        assert self.service._validate_pmid(None) is False

    def test_clean_text(self):
        """Test text cleaning functionality."""
        dirty_text = "This &amp; that &lt;tag&gt; &quot;quoted&quot;"
        clean_text = self.service._clean_text(dirty_text)
        assert clean_text == 'This & that <tag> "quoted"'

        messy_text = "  Multiple   spaces\n\nand\t\ttabs  "
        clean_text = self.service._clean_text(messy_text)
        assert clean_text == "Multiple spaces and tabs"

    def test_get_statistics(self):
        """Test statistics collection."""
        stats = self.service.get_statistics()

        required_keys = [
            "total_searches",
            "total_papers_fetched",
            "total_api_calls",
            "average_response_time",
            "error_count",
            "rate_limit_hits",
        ]

        for key in required_keys:
            assert key in stats
            assert isinstance(stats[key], (int, float))

    def test_reset_statistics(self):
        """Test statistics reset."""
        with patch.object(self.service, "search_papers") as mock_request:
            mock_request.return_value = []
            self.service.search_papers("test")

        self.service.reset_statistics()
        stats = self.service.get_statistics()

        assert stats["total_searches"] == 0
        assert stats["total_api_calls"] == 0

    @patch("src.pubmed_miner.services.paper_collection.Entrez.read")
    @patch("src.pubmed_miner.services.paper_collection.Entrez.esearch")
    def test_search_with_filters(self, mock_esearch, mock_read):
        """Test search with additional filters."""

        def side_effect(handle):
            return _parse_esearch_xml(handle.read())

        mock_read.side_effect = side_effect
        mock_esearch.return_value = _make_esearch_handle([])

        self.service.search_papers(
            "machine learning", date_from="2020/01/01", date_to="2023/12/31"
        )

        call_args = mock_esearch.call_args[1]
        assert "2020/01/01" in call_args["term"]
        assert "2023/12/31" in call_args["term"]

    def test_concurrent_requests(self):
        """Test handling of concurrent requests."""
        import threading

        results = []
        errors = []

        def make_request():
            try:
                with patch.object(self.service, "search_papers") as mock_request:
                    mock_request.return_value = ["12345"]
                    result = self.service.search_papers("test")
                    results.append(result)
            except Exception as e:
                errors.append(e)

        threads = []
        for i in range(5):
            thread = threading.Thread(target=make_request)
            threads.append(thread)
            thread.start()

        for thread in threads:
            thread.join()

        assert len(errors) == 0
        assert len(results) == 5
