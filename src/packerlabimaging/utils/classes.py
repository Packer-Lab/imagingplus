from typing import TypedDict


class ObjectClassError(Exception):
    """handles exceptions caused by calling function on invalid class type."""

    def __init__(self, function, valid_class, invalid_class):
        super().__init__(f'Invalid class ({invalid_class}) being used. <{function}> only available for {valid_class}.')


class IncompatibleFunctionError(Exception):
    """handles exceptions caused by incompatible functions"""

    def __init__(self, function):
        super().__init__(f'Incompatible function call. <{function}>')


class UnavailableOptionError(Exception):
    """handles exceptions caused by unavailable options"""

    def __init__(self, option):
        super().__init__(f'Unavailable option for object. <{option}>')


class PaqInfoTrial(TypedDict, total=False):
    """dictionary with preset keys to hold meta-information about an individual paq file."""

    frame_channel: str
    paq_path: str  # complete path to .paq file for trial
    stim_channel: str


class TrialsInformation(TypedDict, total=False):
    """dictionary with preset keys to hold meta-information about an individual trial."""

    trialType: str
    tiff_path: str
    expGroup: str
    PaqInfoTrial: PaqInfoTrial
    s2p_use: bool
    naparm_path: str
    analysis_object_information: TypedDict("analysis_object_information",
                                           {'series ID': str, 'repr': str, 'pkl path': str})


