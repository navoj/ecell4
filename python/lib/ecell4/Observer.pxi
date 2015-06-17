cdef class Observer:
    """A wrapper for a base class of Observers.

    Warning: This is mainly for developers.
    Do not use this for your simulation.
    """

    def __cinit__(self):
        self.thisptr = new shared_ptr[Cpp_Observer](
            <Cpp_Observer*>(new Cpp_FixedIntervalNumberObserver(
                0.0, vector[string]()))) #XXX: DUMMY

    def __dealloc__(self):
        del self.thisptr

    def next_time(self):
        """Return the next time for logging."""
        return self.thisptr.get().next_time()

    def reset(self):
        """Reset the internal state."""
        self.thisptr.get().reset()

cdef class FixedIntervalNumberObserver:
    """An ``Observer``class to log the number of molecules with the fixed
    step interval.
    This ``Observer`` logs at the current time first, and then keeps logging
    every after the interval.

    FixedIntervalNumberObserver(dt, species)

    """

    def __init__(self, Real dt, species):
        """Constructor.

        Args:
            dt (float): A step interval for logging.
            species (list): A list of strings, but not of ``Species``.
              The strings suggest serials of ``Species`` to be observed.

        """
        pass  # XXX: Only used for doc string

    def __cinit__(self, Real dt, species):
        cdef vector[string] cpp_species
        for serial in species:
            cpp_species.push_back(tostring(serial))
        self.thisptr = new shared_ptr[Cpp_FixedIntervalNumberObserver](
            new Cpp_FixedIntervalNumberObserver(dt, cpp_species))

    def __dealloc__(self):
        del self.thisptr

    def next_time(self):
        """Return the next time for logging."""
        return self.thisptr.get().next_time()

    def num_steps(self):
        """Return the number of steps."""
        return self.thisptr.get().num_steps()

    def data(self):
        """Return a list of the number of molecules you specified.

        Returns:
          list: A list of lists of the numbers of molecules.
            The size of a return value is equal to ``num_steps``.
            Each element of a return value is a list consisting of
            time and the number of molecules specified at the construction.

        """
        cdef vector[vector[Real]] d = self.thisptr.get().data()
        retval = []
        cdef vector[vector[Real]].iterator it = d.begin()
        while it != d.end():
            retval.append(deref(it))
            inc(it)
        return retval

    def targets(self):
        """Return a list of ``Species``, which this ``Observer`` observes

        Returns:
          list: A list of ``Species``. This is generated from arguments
            you gave at the construction.

        """
        cdef vector[Cpp_Species] species = self.thisptr.get().targets()

        retval = []
        cdef vector[Cpp_Species].iterator it = species.begin()
        while it != species.end():
            retval.append(
                 Species_from_Cpp_Species(
                     <Cpp_Species*>(address(deref(it)))))
            inc(it)
        return retval

    def as_base(self):
        """Clone self as a base class. This function is for developers."""
        retval = Observer()
        del retval.thisptr
        retval.thisptr = new shared_ptr[Cpp_Observer](
            <shared_ptr[Cpp_Observer]>deref(self.thisptr))
        return retval

    def reset(self):
        """Reset the internal state."""
        self.thisptr.get().reset()

cdef class NumberObserver:
    """An ``Observer``class to log the number of molecules.
    This ``Observer`` logs at the current time first, and then keeps logging
    every after simulation steps.
    Warning: This doesn't work with ODESimulator.

    NumberObserver(species)

    """

    def __init__(self, species):
        """Constructor.

        Args:
            species (list): A list of strings, but not of ``Species``.
              The strings suggest serials of ``Species`` to be observed.

        """
        pass  # XXX: Only used for doc string

    def __cinit__(self, species):
        cdef vector[string] cpp_species
        for serial in species:
            cpp_species.push_back(tostring(serial))
        self.thisptr = new shared_ptr[Cpp_NumberObserver](
            new Cpp_NumberObserver(cpp_species))

    def __dealloc__(self):
        del self.thisptr

    def num_steps(self):
        """Return the number of steps."""
        return self.thisptr.get().num_steps()

    def next_time(self):
        """Return the next time for logging."""
        return self.thisptr.get().next_time()

    def data(self):
        """Return a list of the numbers of molecules you specified.

        Returns:
          list: A list of lists of the number of molecules.
            The size of a return value is equal to ``num_steps``.
            Each element of a return value is a list consisting of
            time and the number of molecules specified at the construction.

        """
        cdef vector[vector[Real]] d = self.thisptr.get().data()
        retval = []
        cdef vector[vector[Real]].iterator it = d.begin()
        while it != d.end():
            retval.append(deref(it))
            inc(it)
        return retval

    def targets(self):
        """Return a list of ``Species``, which this ``Observer`` observes

        Returns:
          list: A list of ``Species``. This is generated from arguments
            you gave at the construction.

        """
        cdef vector[Cpp_Species] species = self.thisptr.get().targets()

        retval = []
        cdef vector[Cpp_Species].iterator it = species.begin()
        while it != species.end():
            retval.append(
                 Species_from_Cpp_Species(
                     <Cpp_Species*>(address(deref(it)))))
            inc(it)
        return retval

    def as_base(self):
        """Clone self as a base class. This function is for developers."""
        retval = Observer()
        del retval.thisptr
        retval.thisptr = new shared_ptr[Cpp_Observer](
            <shared_ptr[Cpp_Observer]>deref(self.thisptr))
        return retval

    def reset(self):
        """Reset the internal state."""
        self.thisptr.get().reset()

cdef class TimingNumberObserver:
    """An ``Observer``class to log the number of molecules just at the time
    you assigned.

    TimingNumberObserver(t, species)

    """

    def __init__(self, vector[double] t, species):  #XXX: vector[Real]
        """Constructor.

        Args:
            t (list): A list of times for logging. A time prior to the current
              time will be ignored.
            species (list): A list of strings, but not of ``Species``.
              The strings suggest serials of ``Species`` to be observed.

        """
        pass  # XXX: Only used for doc string

    def __cinit__(self, vector[double] t, species):  #XXX: vector[Real]
        cdef vector[string] cpp_species
        for serial in species:
            cpp_species.push_back(tostring(serial))
        self.thisptr = new shared_ptr[Cpp_TimingNumberObserver](
            new Cpp_TimingNumberObserver(t, cpp_species))

    def __dealloc__(self):
        del self.thisptr

    def next_time(self):
        """Return the next time for logging."""
        return self.thisptr.get().next_time()

    def num_steps(self):
        """Return the number of steps."""
        return self.thisptr.get().num_steps()

    def data(self):
        """Return a list of the numbers of molecules you specified.

        Returns:
          list: A list of lists of the number of molecules.
            The size of a return value is equal to ``num_steps``.
            Each element of a return value is a list consisting of
            time and the number of molecules specified at the construction.

        """
        cdef vector[vector[Real]] d = self.thisptr.get().data()
        retval = []
        cdef vector[vector[Real]].iterator it = d.begin()
        while it != d.end():
            retval.append(deref(it))
            inc(it)
        return retval

    def targets(self):
        """Return a list of ``Species``, which this ``Observer`` observes

        Returns:
          list: A list of ``Species``. This is generated from arguments
            you gave at the construction.

        """
        cdef vector[Cpp_Species] species = self.thisptr.get().targets()

        retval = []
        cdef vector[Cpp_Species].iterator it = species.begin()
        while it != species.end():
            retval.append(
                 Species_from_Cpp_Species(
                     <Cpp_Species*>(address(deref(it)))))
            inc(it)
        return retval

    def as_base(self):
        """Clone self as a base class. This function is for developers."""
        retval = Observer()
        del retval.thisptr
        retval.thisptr = new shared_ptr[Cpp_Observer](
            <shared_ptr[Cpp_Observer]>deref(self.thisptr))
        return retval

    def reset(self):
        """Reset the internal state."""
        self.thisptr.get().reset()

cdef class FixedIntervalHDF5Observer:
    """An ``Observer`` class to log the state of ``World`` in HDF5 format
    with the fixed step interval.
    This ``Observer`` saves the ``World`` at the current time first, and
    then keeps saving every after the interval.

    FixedIntervalHDF5Observer(dt, filename)

    """

    def __init__(self, Real dt, filename):
        """Constructor.

        Args:
            dt (float): A step interval for logging.
            filename (str): A file name to be saved. Data are saved in
              HDF5 format. The extension name is recommended to be `.h5`.
              The file name can contain at most one formatting string like
              `%02d`, which will be replaced with the number of steps.
              When the file name contains no formmating string, data will
              be overwritten in a single file at every steps.

        """
        pass  # XXX: Only used for doc string

    def __cinit__(self, Real dt, filename):
        self.thisptr = new shared_ptr[Cpp_FixedIntervalHDF5Observer](
            new Cpp_FixedIntervalHDF5Observer(dt, tostring(filename)))

    def __dealloc__(self):
        del self.thisptr

    def next_time(self):
        """Return the next time for logging."""
        return self.thisptr.get().next_time()

    def num_steps(self):
        """Return the number of steps."""
        return self.thisptr.get().num_steps()

    def filename(self):
        """Return a file name to be saved at the next time"""
        return self.thisptr.get().filename().decode('UTF-8')

    def as_base(self):
        """Clone self as a base class. This function is for developers."""
        retval = Observer()
        del retval.thisptr
        retval.thisptr = new shared_ptr[Cpp_Observer](
            <shared_ptr[Cpp_Observer]>deref(self.thisptr))
        return retval

    def reset(self):
        """Reset the internal state."""
        self.thisptr.get().reset()

cdef class FixedIntervalCSVObserver:
    """An ``Observer`` class to log the state of ``World`` in CSV format
    with the fixed step interval.
    This ``Observer`` saves the ``World`` at the current time first, and
    then keeps saving every after the interval.

    FixedIntervalCSVObserver(dt, filename)

    """

    def __init__(self, Real dt, filename):
        """Constructor.

        Args:
            dt (float): A step interval for logging.
            filename (str): A file name to be saved. Data are saved in
              CSV format. The extension name is recommended to be `.csv`
              or `.txt`.
              The file name can contain at most one formatting string like
              `%02d`, which will be replaced with the number of steps.
              When the file name contains no formmating string, data will
              be overwritten in a single file at every steps.
              The first line in a file represents labels for each row.
              Each column contains a position, a radius, and a serial id
              for the ``Species``.

        """
        pass  # XXX: Only used for doc string

    def __cinit__(self, Real dt, filename, species=None):
        cdef vector[string] cpp_species
        if species is None:
            self.thisptr = new shared_ptr[Cpp_FixedIntervalCSVObserver](
                new Cpp_FixedIntervalCSVObserver(dt, tostring(filename)))
        else:
            for serial in species:
                cpp_species.push_back(tostring(serial))
            self.thisptr = new shared_ptr[Cpp_FixedIntervalCSVObserver](
                new Cpp_FixedIntervalCSVObserver(
                    dt, tostring(filename), cpp_species))

    def __dealloc__(self):
        del self.thisptr

    def next_time(self):
        """Return the next time for logging."""
        return self.thisptr.get().next_time()

    def num_steps(self):
        """Return the number of steps."""
        return self.thisptr.get().num_steps()

    def log(self, w):
        """Force to log the given ``World`` to a file.

        Args:
          w (Space): A ``Space`` (``World``) to be logged.

        Example:
          This is an easy way to save a ``World`` in CSV format without
          running a simulation.

          >>> w = lattice.LatticeWorld(Real3(1, 1, 1), 0.005)
          >>> w.bind_to(NetworkModel())
          >>> w.add_molecules(Species("A"), 3)
          >>> FixedIntervalCSVObserver(1, "test.csv").log(w)
          >>> print(open("test.csv").read())
          x,y,z,r,sid
          0.10614455552060439,0.66106605822212161,0.81500000000000006,0.0050000000000000001,0
          0.38375339303603129,0.37527767497325676,0.23999999999999999,0.0050000000000000001,0
          0.25311394008759508,0.05484827557301445,0.495,0.0050000000000000001,0
        """
        cdef Space space = w.as_base()
        self.thisptr.get().log(space.thisptr.get())

    def filename(self):
        """Return a file name to be saved at the next time"""
        return self.thisptr.get().filename().decode('UTF-8')

    def as_base(self):
        """Clone self as a base class. This function is for developers."""
        retval = Observer()
        del retval.thisptr
        retval.thisptr = new shared_ptr[Cpp_Observer](
            <shared_ptr[Cpp_Observer]>deref(self.thisptr))
        return retval

    def reset(self):
        """Reset the internal state."""
        self.thisptr.get().reset()

cdef class FixedIntervalTrajectoryObserver:
    """An ``Observer`` class to trace and log trajectories of diffusing
    particles in a ``World`` with the fixed step interval.
    This ``Observer`` logs at the current time first, and then keeps logging
    every after the interval.

    FixedIntervalTrajectoryObserver(dt, pids, resolve_boundary=None)

    """

    def __init__(self, Real dt, pids, resolve_boundary=None):
        """Constructor.

        Args:
            dt (float): A step interval for logging.
            pids (list): A list of ``ParticleID``s.
            resolve_boundary (bool, optional): If True, this ``Observer``
              automatically resolves the effect of periodic boundary contidions
              by keeping shifts for each particles. Otherwise, this just
              logs positions within the size of ``World`` with no care about
              boundary conditions.

        """
        pass  # XXX: Only used for doc string

    def __cinit__(self, Real dt, pids, resolve_boundary=None):
        cdef vector[Cpp_ParticleID] tmp
        for pid in pids:
            tmp.push_back(deref((<ParticleID>pid).thisptr))
        if resolve_boundary is None:
            self.thisptr = new shared_ptr[Cpp_FixedIntervalTrajectoryObserver](
                new Cpp_FixedIntervalTrajectoryObserver(dt, tmp))
        else:
            self.thisptr = new shared_ptr[Cpp_FixedIntervalTrajectoryObserver](
                new Cpp_FixedIntervalTrajectoryObserver(dt, tmp, <bool>resolve_boundary))

    def __dealloc__(self):
        del self.thisptr

    def next_time(self):
        """Return the next time for logging."""
        return self.thisptr.get().next_time()

    def num_steps(self):
        """Return the number of steps."""
        return self.thisptr.get().num_steps()

    def data(self):
        """Return a list of trajectories for each particles.

        Returns:
          list: A list of lists of ``Real3``. An element of a return value
            is corresponding the trajectory of each particle. Thus, the size
            of a return value is the same with that of ``pids`` you gave
            at the construction.
            If a particle corresponding to the given ``ParticleID`` is missing,
            i.e. for a reaction, this ``Observer`` just skips to log the
            position. Therefore, lengths of the trajectories can be diverse.

        """
        cdef vector[vector[Cpp_Real3]] d = self.thisptr.get().data()
        retval = []
        cdef vector[vector[Cpp_Real3]].iterator it = d.begin()
        cdef vector[Cpp_Real3].iterator it2
        while it != d.end():
            it2 = deref(it).begin()
            retval.append([])
            while it2 != deref(it).end():
                retval[-1].append(Real3_from_Cpp_Real3(address(deref(it2))))
                inc(it2)
            inc(it)
        return retval

    def as_base(self):
        """Clone self as a base class. This function is for developers."""
        retval = Observer()
        del retval.thisptr
        retval.thisptr = new shared_ptr[Cpp_Observer](
            <shared_ptr[Cpp_Observer]>deref(self.thisptr))
        return retval

    def reset(self):
        """Reset the internal state."""
        self.thisptr.get().reset()
