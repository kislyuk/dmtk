__doc__=""" """
import os
from dmtk.io.SecondaryConfig import SecondaryConfig

SEYMOUR_HOME = None
if 'SEYMOUR_HOME' in os.environ:
    SEYMOUR_HOME = os.environ['SEYMOUR_HOME']

DEFAULT_SMRTANALYSIS_VERSION = '0.0.0'
DEFAULT_SMRTPIPE_VERSION = '0.0.0'

def getSmrtAnalysisVersionInfo( withBuildIds=True, configFilePath=None ):
    if configFilePath is None:
        if not SEYMOUR_HOME or not os.path.exists(SEYMOUR_HOME):
            return DEFAULT_SMRTANALYSIS_VERSION, DEFAULT_SMRTPIPE_VERSION
        configFile = os.path.join( SEYMOUR_HOME, 'etc', 'config.xml' )
    else:
        configFile = configFilePath
    if not os.path.exists(configFile):
        return DEFAULT_SMRTANALYSIS_VERSION, DEFAULT_SMRTPIPE_VERSION
    config = SecondaryConfig( configFile )
    if withBuildIds:
        v1 = '%s.%s' % ( config['version'], config['build'] )
        v2 = '%s.%s' % ( config['SMRTpipe'].version, config['SMRTpipe'].build )
        return v1, v2
    else:
        return config['version'], config['SMRTpipe'].version
