# !usr/bin/env python
# -*- coding: utf-8 -*-
#
# Licensed under a 3-clause BSD license.
#
# @Author: Brian Cherinka
# @Date:   2017-12-05 12:01:21
# @Last modified by:   Brian Cherinka
# @Last Modified time: 2017-12-05 12:19:32

from __future__ import print_function, division, absolute_import


class MugatuError(Exception):
    """A custom core Mugatu exception"""

    def __init__(self, message=None):

        message = 'There has been an error' \
            if not message else message

        super(MugatuError, self).__init__(message)


class MugatuNotImplemented(MugatuError):
    """A custom exception for not yet implemented features."""

    def __init__(self, message=None):

        message = 'This feature is not implemented yet.' \
            if not message else message

        super(MugatuNotImplemented, self).__init__(message)


class MugatuDesignError(MugatuError):
    """A custom exception when there is a critial error in a design."""

    def __init__(self, message=None):

        message = 'There is a critical error in the design.' \
            if not message else message

        super(MugatuDesignError, self).__init__(message)


class MugatuAPIError(MugatuError):
    """A custom exception for API errors"""

    def __init__(self, message=None):
        if not message:
            message = 'Error with Http Response from Mugatu API'
        else:
            message = 'Http response error from Mugatu API. {0}'.format(message)

        super(MugatuAPIError, self).__init__(message)


class MugatuApiAuthError(MugatuAPIError):
    """A custom exception for API authentication errors"""
    pass


class MugatuMissingDependency(MugatuError):
    """A custom exception for missing dependencies."""
    pass


class MugatuWarning(Warning):
    """Base warning for Mugatu."""


class MugatuUserWarning(UserWarning, MugatuWarning):
    """The primary warning class."""
    pass


class MugatuSkippedTestWarning(MugatuUserWarning):
    """A warning for when a test is skipped."""
    pass


class MugatuDesignWarning(MugatuUserWarning):
    """A warning for when a test is skipped."""
    pass


class MugatuDesignModeWarning(MugatuUserWarning):
    """A warning for when a test is skipped."""
    pass


class MugatuDeprecationWarning(MugatuUserWarning):
    """A warning for deprecated features."""
    pass
