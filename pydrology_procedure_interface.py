from .descriptors.list_descriptor import ListDescriptor

class PydrologyProcedureInterface:
    """Abstract class to use as base for pydrology procedures. 
    
    It checks pars, Boundaries and InitialConditions and saves them as object properties. Arguments must be lists, other iterables (i.e. tuple, numpy.array) or None. 
    
    If any argument is not None and not iterable, it will raise a TypeError
    
    Example:

    class MyActualProcedure(PydrologyProcedureInterface):        
        def __init__(pars : Tuple[float, float], Boundaries, InitialConditions):
            super().__init__(pars, Boundaries, InitialConditions)
            ...
    
    Add self.checkItemsType() calls to assert length and type of argument items:

            self.checkItemsType("pars", float, length=2, coerce=True)
    """
    
    pars = ListDescriptor()
    """List of procedure parameters"""

    boundaries = ListDescriptor()
    """List of procedure boundary conditions"""

    initial_conditions = ListDescriptor()
    """List of procedure initial conditions"""

    def __init__(
        self,
        pars : list = None,
        Boundaries : list = None,
        InitialConditions : list = None,
        **kwargs):
        """
        Args:
            pars (List[float]): list of procedure parameters. Defaults to None.
            Boundaries (List[float], optional): list of procedure boundaries. Defaults to None.
            InitialConditions (List[float], optional): list of procedure initial conditions. Defaults to None.
        
        Raises:
            TypeError: if any argument is not None and not iterable
        """

        self.pars = pars
        self.boundaries = Boundaries
        self.initial_conditions = InitialConditions

    def checkItemsType(
        self,
        attr : str,
        item_type : type,
        length : int = None,
        coerce : bool = False
    ) -> None:
        """Iterates attr attribute and raises TypeError if an item is not an instance of item_type. If length is not set, attr may be None

        Args:
            attr (str): attribute name
            item_type (type): type to check items against
            length (int, None): check that length of attr equals this int
            coerce (bool, False): try to force items type and update attr
        """
        try:
            attr_value = getattr(self, attr)
        except AttributeError as e:
            raise AttributeError(e)
        if attr_value is None:
            if length is not None:
                raise TypeError("Attribute %s must be of type %s, not None" % (attr, type))
            else:
                return
        if length is not None and len(attr_value) != length:
            raise ValueError("Attribute %s must be of length %i" % (attr, length))
        for i, v in enumerate(attr_value):
            if not isinstance(v, item_type):
                if coerce:
                    try:
                        new_v = item_type(v)
                    except ValueError as e:
                        raise TypeError("Item %i of %s is not of type %s and can't be coerced" % (i, attr, item_type))
                else:
                    raise TypeError("Item %i of %s is not of type %s" % (i, attr, item_type))
                attr_value[i] = new_v
                