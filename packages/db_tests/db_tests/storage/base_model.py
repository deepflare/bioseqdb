from __future__ import annotations
from sqlalchemy.ext.declarative import declarative_base, DeclarativeMeta
from db_tests.pretty_print import pretty_print_to_str
from pydantic import BaseModel as BaseRawModel, Field
from typing import Iterable, Any, Optional, Type, Dict, List
from db_tests.storage.models_utilities.classproperty import classproperty
from db_tests.storage.fields import OrmFieldsHelper
from db_tests.storage.models_utilities.snake_case import camel_to_snake

def generate_table_name(class_name: str) -> str:
    return f"dw_{camel_to_snake(class_name.replace('Orm', ''))}_model"

def format_str_model(
    description: str,
    header_prefix: str,
    header_title: str,
    values: Dict[str, str],
    values_prefixes: Optional[Dict[str, str]] = None,
    diff_values: Optional[Dict[str, str]] = None,
    no_header: bool = False,
):
    if diff_values:
        diff_values_dict = dict(diff_values)
    global_ident = " "*3
    inner_ident = " "*3
    result: List[str] = [
        "[bold]About:[/bold]\n",
        "|\n",
        *[f"|  {line.strip()}\n" for line in description.split("\n") if len(line.strip()) > 0],
        "|\n",
    ]
    if no_header:
        result = []
    longest_prefix_len = 0
    if values_prefixes:
        longest_prefix_len = max([len(item) for item in values_prefixes.values()])
    for name, value_str in values:
        val_str_part = f"{name}: {value_str}"
        if diff_values:
            val_str_part = f"{name}:\n{global_ident}    [bold red]---[/bold red] {value_str}\n{global_ident}    [bold green]+++[/bold green] {diff_values_dict[name]}"
        if values_prefixes:
            prefix = values_prefixes[name]
            prefix_fill = " " * max(0, longest_prefix_len-len(prefix))
            prefix = prefix.replace("[", "\\[")
            result.append(f"[italic blue]{prefix}[/italic blue]{prefix_fill}   {val_str_part}\n")
        else:
            result.append(f"{val_str_part}\n")
    inner_text = "".join([f"{global_ident}{inner_ident}{line}" for line in result])
    head_text = f"{global_ident}[bold]{header_prefix}[/] [bold red]{header_title}[/]:"
    return f"\n{pretty_print_to_str(head_text)}\n{pretty_print_to_str(inner_text)}"

# Default base model (Pydantic)
class BaseModel(BaseRawModel):
    pass

from typing import Container, Optional, Type

from pydantic import create_model, validator
from sqlalchemy.inspection import inspect
from sqlalchemy.orm.properties import ColumnProperty

IMPORTED_PROPS = [
    *[prop for prop in dir(BaseRawModel) if not prop.startswith("__")],
    "__fields__",
]

def sqlalchemy_to_pydantic(
    db_model: Type, *, exclude: Container[str] = []
) -> Type[BaseModel]:
    class Config:
        orm_mode = True
        # make our setters perform validations (& type coercion) on inputs.
        validate_assignment = True

    mapper = inspect(db_model)
    fields = {}
    validators: Dict[str, Any] = dict()
    for attr in mapper.attrs:
        if isinstance(attr, ColumnProperty):
            if attr.columns:
                name = attr.key
                if name in exclude:
                    continue
                column = attr.columns[0]
                python_type: Optional[type] = None
                if hasattr(column.type, "impl"):
                    if hasattr(column.type.impl, "python_type"):
                        python_type = column.type.impl.python_type
                elif hasattr(column.type, "python_type"):
                    python_type = column.type.python_type
                assert python_type, f"Could not infer python_type for {column}"
                default = None
                if column.default is None and not column.nullable:
                    default = ...
                meta = getattr(column, "_metadata")
                if "validator" in meta:
                    fn = meta["validator"]
                    validators[f"{name}_validator"] = validator(name, pre=True, allow_reuse=True)(fn)
                    del meta["validator"]
                if "validator_post" in meta:
                    fn = meta["validator_post"]
                    def _post_validator(value, values, config, field):
                        if not fn(value, values, config, field):
                            raise ValueError("Invalid value was provided: {value}")
                        return value
                    validators[f"{name}_validator_post"] = validator(name, pre=False, allow_reuse=True)(_post_validator)
                    del meta["validator_post"]
                if "type" in meta:
                    python_type = meta["type"]
                    del meta["type"]
                fields[name] = (python_type, Field(default, **meta))
    pydantic_model = create_model(
        db_model.__name__, __config__=Config, __validators__=validators, **fields  # type: ignore
    )
    return pydantic_model

# Metaclass for SQLAlchemy base model
class OrmBaseMetaModel(DeclarativeMeta):
    def __init__(cls, *args):
        if cls.__name__ != "Base":
            table_name = generate_table_name(cls.__name__)
            cls.__tablename__ = table_name
            cls.__table_args__ = dict()
        super().__init__(*args)
        if cls.__name__ == "Base":
            return

        def _model_validator(v: Any):
            assert isinstance(v, cls)
            return cls(**v.dict())

        print(f"Register {cls.__name__}")
        from pydantic.validators import _VALIDATORS
        _VALIDATORS.append((cls, [_model_validator]))

        if not cls.__doc__ or len(cls.__doc__) == 0:
            raise ValueError(f"Please provide a docstring for a class {cls.__name__}(): We care about the documentation of the models.")
        OrmBaseModelPrototype._OrmBaseModelPrototype__orm_models[cls.__name__] = cls
    
    def __call__(cls, *args, **kwargs):
        obj = type.__call__(cls, *args, **kwargs, _will_have_post_init = True)
        # Post __init__ hook
        obj._post_init(cls, *args, **kwargs)
        return obj

    def __getattr__(cls, name):
        if name in cls.__dict__:
            return cls.__dict__[name]
        if name in IMPORTED_PROPS:
            m = cls.model
            if hasattr(m, name):
                return getattr(m, name)
        raise AttributeError

    def __str__(cls) -> str:
        vals = cls.model.__fields__.items()
        return format_str_model(
            description=cls.__doc__,
            header_prefix="class",
            header_title=cls.__name__,
            values=[(name, f"[green]{value.field_info.description}[/green]") for name, value in vals],
            values_prefixes={name: ("None" if value is None else value._type_display()) for name, value in vals},
        )

# SQLAlchemy base model
class OrmBaseModelPrototype:
    __tablename__ = generate_table_name("EmptyOrm")
    __table_args__ = {'extend_existing': True}

    __orm_models: Dict[str, Type[OrmBaseModel]] = dict()
    _pydantic_model: Optional[Type[BaseModel]] = None
    _model_values: Optional[BaseModel] = None
    _will_have_post_init: bool = False
    _frozen: bool = False

    # @classproperty
    # def query(cls):
    #     return STORAGE.db.query(cls)

    @classproperty
    def model(cls) -> Type[BaseModel]:
        if cls._pydantic_model is None:
            cls._pydantic_model = sqlalchemy_to_pydantic(cls)
        return cls._pydantic_model

    @classproperty
    def __pydantic_model__(cls) -> Type[BaseModel]:
        return cls.model

    @property
    def model_values(self) -> BaseModel:
        if self._model_values is None:
            self._model_values = self.model.from_orm(self)
        return self._model_values

    def _post_init(self, cls: Type[OrmBaseModel], *args, **kwargs):
        self._model_values = self.model.from_orm(self)

    # @classmethod
    # def bulk_create(
    #     cls,
    #     values: Iterable[OrmBaseModel],
    # ):
    #     STORAGE.db.execute(
    #         cls.__table__.insert(),
    #         [item.model_values.dict() for item in values],
    #     )
    
    @classmethod
    def get_comparator(
        cls,
        ignore_attrs: Optional[List[str]] = None,
    ):
        def _cmp(x, y):
            if x.__class__.__name__ != y.__class__.__name__:
                return False
            for name, val in x.model_values.__repr_args__():
                if name in ignore_attrs:
                    # Ignore
                    continue
                eq = getattr(x.model_values, name) == getattr(y.model_values, name)
                if not eq:
                    return False
            return True
        return _cmp

    def __eq__(self, other):
        x, y = self, other
        
        ignore_attrs = []
        if hasattr(x.__class__, "__ignore_attrs_eq__"):
            ignore_attrs += getattr(x.__class__, "__ignore_attrs_eq__")
        if hasattr(y.__class__, "__ignore_attrs_eq__"):
            ignore_attrs += getattr(y.__class__, "__ignore_attrs_eq__")
        ignore_attrs = set(ignore_attrs)
        
        if x.__class__.__name__ != y.__class__.__name__:
            return False
        for name, val in x.model_values.__repr_args__():
            if name in ignore_attrs:
                # Ignore
                continue
            eq = getattr(x.model_values, name) == getattr(y.model_values, name)
            if not eq:
                return False
        return True
        #if hasattr(other, "model_values"):
        #    return self.model_values.__eq__(other.model_values)
        #return False

    def to_str(
        self,
        diff_from = None,
        no_header: bool = False,
    ) -> str:
        vals = self.model_values
        model_fields = self.__class__.model.__fields__.items()

        values = []
        diff_values = None
        if diff_from:
            diff_values = []
            other_vals_dict = dict(diff_from.model_values.__repr_args__())
            for name, value in vals.__repr_args__():
                other_value = other_vals_dict[name] 
                if value != other_value:
                    values.append((name, ("None" if value is None else repr(value))))
                    diff_values.append((name, ("None" if other_value is None else repr(other_value))))
        else:
            for name, value in vals.__repr_args__():
                values.append((name, ("None" if value is None else repr(value))))

        values_prefixes={name: ("None" if value is None else f"{value._type_display()}") for name, value in model_fields}
        return format_str_model(
            description=self.__class__.__doc__,
            header_prefix="object",
            header_title=self.__class__.__name__,
            values=values,
            diff_values=diff_values,
            values_prefixes=values_prefixes,
            no_header=no_header,
        )

    def __str__(self) -> str:
        return self.to_str()

    def json(self, *args, **kwargs):
        return self.model_values.json(*args, **kwargs)

def _custom_model_getattr(self, name):
    if name in self.__dict__:
        return self.__dict__[name]
    if name in IMPORTED_PROPS:
        m = self.model_values
        if hasattr(m, name):
            return getattr(m, name)
    raise AttributeError

def _custom_model_setattr(self, name, value):
    if name == "_model_values" and self._frozen:
        pass
    elif name not in ["_will_have_post_init", "_frozen"] and (self._model_values is not None or self._frozen):
        raise Exception(f"Cannot modify attribute on class OrmBaseModel: {name} => {value}")
    else:
        #print(f"Allow change {name} => {value}")
        if name == "_will_have_post_init" and value == True:
            self._frozen = False
    self.__dict__[name] = value
    if name == "_sa_instance_state" and not self._will_have_post_init:
        self._frozen = True

# def _custom_model_create_all():
#      STORAGE._init_db(OrmBaseModel)

OrmBaseModel = declarative_base(cls=OrmBaseModelPrototype, metaclass=OrmBaseMetaModel)
org_init = OrmBaseModel.__init__

OrmBaseModel.__setattr__ = _custom_model_setattr
OrmBaseModel.__getattr__ = _custom_model_getattr
#OrmBaseModel.create_all = _custom_model_create_all

def lmao(self, *args, _will_have_post_init = False, **kwargs):
    self._will_have_post_init = True
    return org_init(self, *args, **kwargs)

OrmBaseModel.__init__ = lmao

FIELDS = OrmFieldsHelper()