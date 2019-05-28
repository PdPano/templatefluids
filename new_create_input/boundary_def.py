"""
Defines boundary container
"""

class Boundary:
    def __init__(self, boundary_dict):
        # Mandatory keys
        try:
            self.type = boundary_dict['type']
            self.position = boundary_dict['position']
            self.xlims = boundary_dict['xlims']
            self.ylims = boundary_dict['ylims']
        except KeyError as exc:
            print("***************")
            print("Key {} is mandatory in BOUNDARY entry".format(exc.args[0]))
            print("***************")
            raise
        # Optional parameters
        self.timeFunctionType = boundary_dict.get("timeFunctionType", "0.0")
        self.timeParameterA = boundary_dict.get("timeParameterA", "1.0")
        self.timeParameterB = boundary_dict.get("timeParameterB", "0.0")
