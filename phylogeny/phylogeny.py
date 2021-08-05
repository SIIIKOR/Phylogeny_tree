"""
Mateusz Sikorski

"""
from Bio import SeqIO, Seq
from numpy import zeros
from collections import deque
from random import choice

import datetime
from typing import List, Tuple


class Sample:
    sample_id: str
    country: str
    date: datetime.date
    seq: Seq.Seq    # Jeśli ktoś woli, może być też Bio.Seq.Seq z Biopythona.

    def __init__(self, description=None, seq=None):
        """
        wierzchołek drzewa

        posiada liste dzieci
        :param str description: zawiera sample_id, counry, date
        :param Bio.Seq.Seq seq: sekwencja
        """
        if description is not None and seq is not None:  # żeby dało się tworzyć obiekty z settera
            self.sample_id, self.country, self.date = description.split("|")
            self.date = datetime.datetime.strptime(self.date, "%Y-%m-%d").date()
            self.seq = seq
            self.children = []
        self.parent = None  # atrybut wskazujący na ojca, użyteczne tylko w funckji construct_approximate_tree

    def set_parent(self, parent):
        """
        ustawia ojca dla tego nodea, tak naprawdę metoda używanan tylko w construct_approximate_tree

        :param sample parent:
        :return: None
        """
        self.parent = parent

    def add_child(self, sample):
        """
        dodaje dziecko do tego wierzchołka

        :param Sample sample:
        :return: None
        """
        self.children.append(sample)

    def set_attributes(self, sample_id, country, date, seq, parent=None):
        """
        setter wykorzystany do kopiowania obiektów sample, wiem, że jest coś takiego jak deepcopy()
        ale nie wiem czy jest mi potrzebny

        :param str sample_id:
        :param str country:
        :param datetime.date date:
        :param Seq.Seq seq:
        :param sample parent:
        :return: None
        """
        self.sample_id = sample_id
        self.country = country
        self.date = date
        self.seq = seq
        self.parent = parent
        self.children = []

    def get_copy(self):
        """
        :return: kopia self, która domyślnie nie ma ustawionych dzieci
        """

        new_sample = Sample()
        new_sample.set_attributes(self.sample_id, self.country, self.date, self.seq, self.parent)
        return new_sample

    def neighbours(self):
        """
        metoda potrzebna do construct_approximate_tree
        :return: sąsiedzi nodea
        """
        if self.parent is None:
            return self.children
        return self.children + [self.parent]  # pojedyńcze dodaje do większego - lab6


class Tree:
    def __init__(self, root, cost=None):
        self.root = root
        self.cost = cost

    def edges(self) -> List[Tuple[str, str]]:
        """
        zastosowałem algorytm bfs, ale dfs jest równie dobry

        drobna uwaga, skoro to są drzewa nie muszę pamiętać odwiedzonych nodeów,
        bo z definicji każdy odwiedzam raz
        :return: edges jak w definicji
        """
        # TODO: funkcja powinna zwracać listę wszystkich krawędzi w drzewie,
        # w dowolnej kolejności. Każda krawędź to para (id_próbki_rodzica, id_próbki_dziecka).
        # Wymagana złożoność: O(n), gdzie n to liczba wierzchołków w drzewie.
        queue = deque([self.root])
        edges = []
        while queue:
            node = queue.pop()
            for child in node.children:
                queue.appendleft(child)
                edges.append((node.sample_id, child.sample_id))
        return edges

    def filter(self, country: str) -> List['Tree']:
        """
        funkcja stworzona na bazie dfs

        :param str country:
        :return: lista Tree
        """
        # TODO: funkcja powinna zwracać listę nowych drzew,
        # składających się z próbek pochodzących z podanego kraju,
        # według definicji z treści zadania.
        # Wymagana złożoność: O(n), gdzie n to liczba wierzchołków w drzewie.
        stack = [(self.root, None)]  # każdy node musi pamiętać ostatniego pasującego ojca, jęsli go nie ma to None
        trees = []
        while stack:
            node, parent = stack.pop()
            if node.country == country:
                copied_node = node.get_copy()  # musimy stworzyć kopię, żeby nie było problemów z dziećmi
                if parent is None:  # jeśli zaczynamy nowe poddrzewo
                    trees.append(Tree(copied_node))  # tworzymy drzewo i node jest korzeniem
                else:
                    parent.add_child(copied_node)  # przeciwnie dodajemy node do ostatniego pasującego ojca
                parent = copied_node  # ostatni pasujący ojciec
            for child in node.children:
                stack.append((child, parent))  # dziecko, ojciec
        return trees


def read_data(filename: str) -> List[Sample]:
    """
    wykorzystałem biopythona

    :param filename: nazwa dostarczonego pliku
    :return: posortowane rosnąco sample w liscie
    """
    # TODO: funkcja powinna wczytywać informacje o próbkach z pliku FASTA o podanej nazwie,
    # tworzyć obiekty typu Sample i zwrócić je, posortowane rosnąco według daty pobrania próbki.
    # Oczekiwana złożoność: O(n * m + n log n), gdzie n to liczba próbek, a m to maksymalna długość sekwencji.
    # Plik FASTA składa się z wielu sekwencji poprzedzonych opisami postaci
    # >ID_próbki|kraj|data
    # (ID próbki, nazwa kraju ani data nie zawierają znaku "|").
    # Na przykład: >MT066156.1|Italy|2020-3-13 (patrz: przykładowe testy).
    # Można wykorzystać bibliotekę Biopython do łatwiejszego wczytania pliku.
    samples = []
    for record in SeqIO.parse(filename, "fasta"):
        samples.append(Sample(record.description, record.seq))
    return sorted(samples, key=lambda sample: sample.date)


def construct_optimal_tree(samples: List[Sample]) -> Tree:
    """
    dla każdego samplea począwszy od drugiego, bo pierwszy jest zawsze rootem,
    szukam optymalnego rodzica. Przy okazji zliczam koszt drzewa

    :param samples: posortowana rosnąco lista sampli
    :return: optymalne drzewo
    """
    # TODO: funkcja powinna zwracać optymalne drzewo filogenetyczne dla podanej listy próbek.
    # Wymagana złożoność: O(n^2 * m^2), gdzie n to liczba próbek, a m to maksymalna długość sekwencji.
    total_score = 0
    for child_index in range(1, len(samples)):  # w prawo
        best = None
        for parent_index in range(child_index-1, -1, -1):  # w lewo
            score = levenshtein(samples[child_index].seq, samples[parent_index].seq)  # liczę score
            if best is None or score < best[1]:
                best = (samples[parent_index], score)  # najbardziej pasującu ojciec i score dopasowania
        best[0].add_child(samples[child_index])  # do ojca podłączam dziecko
        total_score += best[1]
    return Tree(samples[0], total_score)


def construct_approximate_tree(samples: List[Sample]) -> Tree:
    """
    :param samples:
    :return: approximated tree
    """
    # TODO: funkcja powinna zwracać drzewo filogenetyczne dla podanej listy próbek.
    # Funkcja powinna działać wyraźnie szybciej niż construct_optimal_tree;
    # wymagany orientacyjny koszt wynikowego drzewa można znaleźć w plikach .txt.
    root = samples[0]  # najstarszy sample jest korzeniem
    all_nodes = [root]  # zbiór wierzchołków w drzewie
    scores = modified_dict()  # score to distance, czym niższy tym lepszy
    total_score = 0
    for child_index in range(1, len(samples)):  # dla każdego nodea poza rootem
        candidates = []
        parent, score = None, None
        for _ in range(3):  # mieliśmy zbierać po 3 kandydatów
            stack = [choice(all_nodes)]  # losowy ojciec ze zbioru wierzchołków w drzewie
            # stack dlatego, bo trzeba obsłużyć potencjalną sytuacje, że jakiś sąsiad ojca, ma mniejszy score
            # z dzieckiem, wtedy trzeba ustawić tego sąsiada jako "nowego" ojca i powtórzyć operacje
            while stack:
                parent = stack.pop()
                score = scores.get_score(samples[child_index], parent)
                neighbour_scores = []  # zbieram wyniki z sąsiadami
                for neighbour in parent.neighbours():
                    new_score = scores.get_score(samples[child_index], neighbour)
                    if new_score < score:
                        neighbour_scores.append((neighbour, new_score))
                if len(neighbour_scores) != 0:
                    stack.append(min(neighbour_scores, key=lambda x: x[1])[0])  # wybieram sąsiada o najmniejszym score
            candidates.append((parent, score))
        best_candidate = min(candidates, key=lambda candidate: candidate[1])  # wybieram najlepszego kandydata
        samples[child_index].set_parent(best_candidate[0])  # ustawiam rodzica dziecka
        best_candidate[0].add_child(samples[child_index])  # dodaje dziecko do rodzica
        all_nodes.append(samples[child_index])  # dodaje dziecko do zbioru wszystkich wierzchołków
        total_score += best_candidate[1]  # powiększam koszt drzewa o dodany wierzchołek

    return Tree(root, total_score)


class modified_dict:
    def __init__(self):
        """
        ten "zmodyfikowany" słownik modelowo jako klucze przyjmuje "zbiory"
        czyli krotki sample1, sample2 i sample2, sample1 powinny odpowiadać
        temu samemu kluczowi i słownik zapytany o score dla obu krotek powinien
        dawać ten sam score

        robię to dlatego, żeby nie liczyc dwa razy odległości edycyjnej

        uwaga działa tylko dla krotek dwuelementowych
        """
        self.normal_dict = {}

    def get_score(self, a, b):
        """
        dostaje score dla dwóch samplei, zapamiętuje wcześniej policzone
        gdy dla jakieś pary nie zostały policzone to je liczy i zwraca
        :param sample a:
        :param sample b:
        :return:
        """
        if (a, b) in self.normal_dict:
            score = self.normal_dict[(a, b)]
        elif (b, a) in self.normal_dict:
            score = self.normal_dict[(b, a)]
        else:
            score = levenshtein(a.seq, b.seq)
            self.normal_dict[(a, b)] = score
        return score


def levenshtein(a, b):
    """
    :param Seq.Seq a:
    :param Seq.Seq b:
    :return: odległość edycyjna
    """
    match = 0
    mis_match = 1
    gap = 1

    m = zeros((len(a)+1, len(b)+1), dtype=int)
    for x in range(1, len(b)+1):  # niestety musze "ręcznie" ustawić pierwszą kolumne i pierwszy wiersz
        m[0, x] = m[0, x-1] + 1
    for y in range(1, len(a)+1):
        m[y, 0] = m[y-1, 0] + 1

    for y in range(1, len(a)+1):
        for x in range(1, len(b)+1):
            left = m[y-1, x] + gap
            right = m[y, x-1] + gap
            if a[y-1] == b[x-1]:
                diag = m[y-1, x-1] + match
            else:
                diag = m[y-1, x-1] + mis_match
            chosen = min(diag, left, right)
            m[y, x] = chosen

    return m[len(a), len(b)]
