import pytest

from koel.cli import main
from koel.sources import SOURCES


def test_download_list_runs(capsys):
    assert main(["download", "--list"]) == 0
    assert "civic" in capsys.readouterr().out


def test_every_source_declares_a_licence():
    assert SOURCES
    for s in SOURCES:
        assert s.licence and s.url


def test_unknown_command_fails():
    with pytest.raises(SystemExit):
        main(["not-a-command"])
