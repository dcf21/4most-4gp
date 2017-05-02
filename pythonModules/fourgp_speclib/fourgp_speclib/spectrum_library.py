#!/usr/bin/env python
# -*- coding: utf-8 -*-


def requires_ids_or_filenames(method):
    """
    A decorator for spectrum library methods that require either a list of Ids or a list of filenames.

    :param method:
        A method belonging to a sub-class of SpectrumLibrary.
    """

    def wrapper(model, *args, **kwargs):
        have_ids = ("ids" in kwargs) and (kwargs["ids"] is not None)
        have_filenames = ("filenames" in kwargs) and (kwargs["filenames"] is not None)
        assert have_ids or have_filenames, "Must supply a list of Ids or a list of filenames"
        assert not (have_ids and have_filenames), "Must supply either a list of Ids or a list of filenames, not both."

        # If a single Id is supplied, rather than a list of Ids, turn it into a one-entry tuple
        if have_ids and not isinstance(kwargs["ids"], (list, tuple)):
            kwargs["ids"] = (kwargs["ids"],)

        # If a single filename is supplied, turn it into a one-entry tuple
        if have_filenames and not isinstance(kwargs["filenames"], (list, tuple)):
            kwargs["filenames"] = (kwargs["filenames"],)

        return method(model, *args, **kwargs)

    return wrapper


class SpectrumLibrary(object):
    """
    An abstract spectrum library object
    """

    def __str__(self):
        return "<{module}.{name} instance".format(module=self.__module__,
                                                  name=type(self).__name__)

    def __repr__(self):
        return "<{0}.{1} object at {2}>".format(self.__module__,
                                                type(self).__name__, hex(id(self)))
