import random
from enum import Enum
import string
from typing import List, Optional

def random_string(
    alphabet_str: str,
    length: int,
) -> str:
    return ''.join(random.SystemRandom().choice(alphabet_str) for _ in range(length))

class Alphabet(str, Enum):
    NUCL_SEQ = "ACTG" #"ACGTX"
    AA_SEQ = "ARNDCQEGHILKMFPSTWYV*"
    AMB_AA_SEQ = "ARNDCQEGHILKMFPSTWYVBJZX*"
    ANY = string.printable

    def random_list(
        self,
        min_length: int,
        max_length: int,
        min_count: int,
        max_count: int,
    ) -> List[str]:
        return [self.random(length=random.randint(min_length, max_length)) for _ in range(random.randint(min_count, max_count))]

    def random(self, length: int) -> str:
        return random_string(
            alphabet_str=self.value,
            length=length,
        )

class SequenceGenerator:
    def __init__(self) -> None:
        pass

    def __call__(
        self,
        description: str,
        input_alphabet: Alphabet,
        min_length: int = 0,
        length_variation: float = 1000.0,
        min_count: int = 10,
        count_variation: float = 100.0,
        max_length: Optional[int] = None,
        max_count: Optional[int] = None,
    ) -> List[str]:
        if max_count is None:
            max_count = min_count + int(count_variation/100.0 * min_count)
        if max_length is None:
            max_length = min_length + int(length_variation/100.0 * min_length)
        return input_alphabet.random_list(
            min_length=min_length,
            max_length=max_length,
            min_count=min_count,
            max_count=max_count,
        )