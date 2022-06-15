import random
from enum import Enum
import string
from typing import List

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