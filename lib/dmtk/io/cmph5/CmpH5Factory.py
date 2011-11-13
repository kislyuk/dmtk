import os
import h5py

##################################################################
#
# This factory is the safest/best way to generate a CmpH5 object.
#
class Singleton(type):
    """Standard Singleton pattern"""
    def __init__( cls, name, bases, dict ):
        super( Singleton, cls ).__init__( name, bases, dict )
        cls.instance = None
    def __call__( cls, *args, **kwargs ):
        if cls.instance == None:
            cls.instance = super( Singleton, cls ).__call__( *args, **kwargs )
        return cls.instance

class CmpH5Factory( object ):
    """This class generates CmpH5 objects (Astro, SF, etc.) from filename/paths."""
    __metaclass__ = Singleton

    def __init__( self ):
        self._opened = set( )

    def destroy( self, cmpH5 ):
        """Signals the factory that the file has been closed and we can allow it to
        be opened again."""
        try:
            self._opened.remove( cmpH5.filename )
        except KeyError:
            pass

    def create( self, cmpH5Path, mode='r', cmpType=None, readType='standard' ):
        """Given the path to a valid .cmp.h5 file, returns a CmpH5 object. 
        If the mode is read or append, the type of CmpH5 to create is determined from
        the version string in the file. If the mode is write, then the optional
        'type' parameter to this function is used.
        This is the best way to create CmpH5 objects."""

        if mode in "ra":
            cmpH5 = h5py.File( cmpH5Path, 'r' )
            version = str(cmpH5.attrs["Version"]) if "Version" in cmpH5.attrs else None
            readType = str(cmpH5.attrs["ReadType"]) if "ReadType" in cmpH5.attrs else None
            from CmpH5 import SFCmpH5
            from CmpH5 import PBCmpH5
            from CmpH5 import SFRCCSCmpH5
            from __init__ import V1_2_0_SF, V1_2_0_PB, V1_2_0
            v2c = { V1_2_0_SF: SFCmpH5,
                    V1_2_0_PB: PBCmpH5,
                    V1_2_0   : SFCmpH5 }

            if readType == "RCCS" and version == V1_2_0:
                cmpType = SFRCCSCmpH5
            else:
                cmpType = SFCmpH5 if version == None else v2c[ version.strip() ]

            cmpH5.close()
        
        if cmpType == None:
            raise IOError, "A CmpH5 Type must be specified when creating a new cmp.h5 file. (%s)" % cmpH5Path
        cmpH5 = cmpType( os.path.abspath( cmpH5Path ), mode=mode, readType=readType )

        if cmpH5.filename in self._opened:
            raise IOError, "An opened CmpH5 object already exists for this file (%s)" % cmpH5.filename

        self._opened.add( cmpH5.filename )
        return cmpH5

factory = CmpH5Factory( )
