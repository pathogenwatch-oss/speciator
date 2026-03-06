import pytest

from speciator.linking import LinkedGroup, Linker


@pytest.fixture
def linked_group():
    return LinkedGroup()

def test_init(linked_group):
    assert len(linked_group.members) == 0
    assert len(linked_group.scores) == 0

def test_len(linked_group):
    assert len(linked_group) == 0
    linked_group.add_score("id1", "id2", 0.5)
    assert len(linked_group) == 2

def test_contains_id(linked_group):
    assert not linked_group.contains_id("id1")
    linked_group.add_score("id1", "id2", 0.5)
    assert linked_group.contains_id("id1")
    assert linked_group.contains_id("id2")

def test_add_score(linked_group):
    linked_group.add_score("id1", "id2", 0.5)
    assert len(linked_group.members) == 2
    assert len(linked_group.scores) == 1
    assert ("id1", "id2") in linked_group.scores
    assert linked_group.scores[("id1", "id2")] == 0.5

    # Test adding score with reversed order of ids
    linked_group.add_score("id3", "id2", 0.7)
    assert len(linked_group.members) == 3
    assert len(linked_group.scores) == 2
    assert ("id2", "id3") in linked_group.scores
    assert linked_group.scores[("id2", "id3")] == 0.7

def test_merge_groups():
    linker = Linker()
    
    # Create two separate groups
    linker.add_score("id1", "id2", 0.5)
    linker.add_score("id3", "id4", 0.7)
    
    assert len(linker.groups) == 2
    
    # Merge the groups
    linker.add_score("id2", "id3", 0.6)
    
    assert len(linker.groups) == 1
    
    merged_group = next(iter(linker.groups))
    
    # Check if all members are in the merged group
    assert merged_group.contains_id("id1")
    assert merged_group.contains_id("id2")
    assert merged_group.contains_id("id3")
    assert merged_group.contains_id("id4")
    
    # Check if all scores are preserved
    assert merged_group.scores[("id1", "id2")] == 0.5
    assert merged_group.scores[("id3", "id4")] == 0.7
    assert merged_group.scores[("id2", "id3")] == 0.6
    
    # Check if all ids are mapped to the merged group
    assert all(linker.id_to_group[name] is merged_group for name in ["id1", "id2", "id3", "id4"])
