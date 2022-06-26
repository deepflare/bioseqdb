from typing import Callable, List, Optional
from db_tests.storage.base_model import OrmBaseModel
import pytest_diff
from itertools import zip_longest

class RecordsList:

    content: List[OrmBaseModel]
    custom_comparator: Optional[Callable[..., bool]]
    ignore_order: bool

    def __init__(
        self,
        content: List[OrmBaseModel],
        compare: Optional[Callable[..., bool]] = None,
        ignore_order: bool = False,
    ) -> None:
        self.content = content
        self.custom_comparator = compare
        self.ignore_order = ignore_order
    
    def __eq__(self, __o: object) -> bool:
        if not isinstance(__o, RecordsList):
            return False
        for index, (x, y) in enumerate(zip_longest(self.content, __o.content)):
            #for index, item in enumerate(self.content):
            eq = False
            if self.custom_comparator:
                eq = self.custom_comparator(x, y)
            else:
                eq = (x == y)
            if not eq:
                if self.ignore_order or __o.ignore_order:
                    eq = False
                    for second_item in __o.content:
                        if self.custom_comparator:
                            eq = self.custom_comparator(x, second_item)
                        else:
                            eq = (x == second_item)
                        if eq:
                            break
                    if not eq:
                        return False
                    else:
                        continue
                return False
        return True
    
    def __repr__(self) -> str:
        return f"RecordsList({len(self.content)} items)"
        # output: List[str] = []
        # output.append(f"RecordsList ({len(self.content)} items) {{")
        # for index, item in enumerate(self.content):
        #     output.append(f"  {index} ==> {item}")
        # output.append("}")
        # return "\n".join(output)

MAX_DIFF_MISMATCH_COUNT = 5

@pytest_diff.registry.register(RecordsList)
def diff(x: RecordsList, y: RecordsList):
    diff_str = []
    mismatches = []
    mismatch_no = 0
    for index, (item_x, item_y) in enumerate(zip_longest(y.content, x.content)):
        if x.ignore_order or y.ignore_order:
            eq = False
            for second_item in y.content:
                eq = (item_x == second_item)
                if eq:
                    break
            if not eq:
                mismatches.append((mismatch_no, index, item_x, item_y))
                mismatch_no = mismatch_no + 1
        else:
            if item_x != item_y:
                mismatches.append((mismatch_no, index, item_x, item_y))
                mismatch_no = mismatch_no + 1
    for mismatch_no, original_pos, x, y in mismatches:
        if x is None or y is None:
            diff_str.append(f"Record #{original_pos} mismatch:")
            if x is not None:
                diff_str = diff_str + x.to_str(no_header=(mismatch_no > 0)).splitlines()
            else:
                diff_str.append("Missing record: None")
            if y is not None:
                diff_str = diff_str + y.to_str(no_header=(mismatch_no > 0)).splitlines()
            else:
                diff_str.append("Missing record: None")
        else:
            diff_str.append(f"Record #{original_pos} mismatch:")
            diff_str = diff_str + y.to_str(diff_from=x, no_header=(mismatch_no > 0)).splitlines()
            if mismatch_no > MAX_DIFF_MISMATCH_COUNT:
                diff_str.append("")
                diff_str.append(f"Rows diffs were limited to {MAX_DIFF_MISMATCH_COUNT} records.")
                diff_str.append(f"There are {len(mismatches)-mismatch_no-1} more different rows.")
                break
    return diff_str
