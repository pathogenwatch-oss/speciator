from __future__ import annotations


class LinkedGroup:
    def __init__(self):
        self.members: set[str | int] = set()
        self.scores: dict[tuple[str | int, str | int], float] = {}

    def __len__(self):
        return len(self.members)

    def contains_id(self, identifier: str | int) -> bool:
        return identifier in self.members

    def add_score(self, id1: str | int, id2: str | int, score: float):
        self.members.add(id1)
        self.members.add(id2)
        self.scores[tuple(sorted([id1, id2]))] = score


class Linker:
    """NB Not thread safe."""

    def __init__(self):
        self.groups: set[LinkedGroup] = set()
        self.id_to_group: dict[str | int, LinkedGroup] = {}

    def add_scores_batch(
        self, scores: list[tuple[str | int, str | int, float]]
    ) -> None:
        for id1, id2, score in scores:
            self.add_score(id1, id2, score)

    def add_score(self, id1: str | int, id2: str | int, score: float) -> None:
        group1 = self.id_to_group.get(id1)
        group2 = self.id_to_group.get(id2)

        if group1 is group2:
            if group1:
                group1.add_score(id1, id2, score)
            else:
                new_group = LinkedGroup()
                new_group.add_score(id1, id2, score)
                self.groups.add(new_group)
                self.id_to_group[id1] = new_group
                self.id_to_group[id2] = new_group
        elif group1 and group2:
            self.merge_groups(group1, group2)
            group1.add_score(id1, id2, score)
        elif group1:
            group1.add_score(id1, id2, score)
            self.id_to_group[id2] = group1
        elif group2:
            group2.add_score(id1, id2, score)
            self.id_to_group[id1] = group2

    def merge_groups(self, group1: LinkedGroup, group2: LinkedGroup):
        if len(group1) < len(group2):
            group1, group2 = group2, group1

        group1.members.update(group2.members)
        group1.scores.update(group2.scores)

        for member in group2.members:
            self.id_to_group[member] = group1

        self.groups.remove(group2)
