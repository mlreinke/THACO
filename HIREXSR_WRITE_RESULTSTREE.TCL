ADD NODE \SPECTROSCOPY::TOP.HIREXSR.RESULTS.BSTI
ADD NODE/USAGE=NUMERIC \SPECTROSCOPY::TOP.HIREXSR.RESULTS.BSTI:CONFIG
ADD NODE/USAGE=SIGNAL \SPECTROSCOPY::TOP.HIREXSR.RESULTS.BSTI:DATA
ADD NODE/USAGE=SIGNAL \SPECTROSCOPY::TOP.HIREXSR.RESULTS.BSTI:FIT
ADD NODE/USAGE=SIGNAL \SPECTROSCOPY::TOP.HIREXSR.RESULTS.BSTI:INST
SET NODE/MODEL_WRITE/SHOT_WRITE/NOWRITE_ONCE \SPECTROSCOPY::TOP.HIREXSR.RESULTS.BSTI/COMPRESS_ON_PUT/NODO_NOT_COMPRESS
SET NODE /ON/NOINCLUDE \SPECTROSCOPY::TOP.HIREXSR.RESULTS.BSTI
SET NODE/MODEL_WRITE/SHOT_WRITE/NOWRITE_ONCE \SPECTROSCOPY::TOP.HIREXSR.RESULTS.BSTI:CONFIG/COMPRESS_ON_PUT/NODO_NOT_COMPRESS
PUT \SPECTROSCOPY::TOP.HIREXSR.RESULTS.BSTI:CONFIG ""
SET NODE /ON/NOINCLUDE \SPECTROSCOPY::TOP.HIREXSR.RESULTS.BSTI:CONFIG
SET NODE/MODEL_WRITE/SHOT_WRITE/NOWRITE_ONCE \SPECTROSCOPY::TOP.HIREXSR.RESULTS.BSTI:DATA/COMPRESS_ON_PUT/NODO_NOT_COMPRESS
PUT \SPECTROSCOPY::TOP.HIREXSR.RESULTS.BSTI:DATA ""
SET NODE /ON/NOINCLUDE \SPECTROSCOPY::TOP.HIREXSR.RESULTS.BSTI:DATA
SET NODE/MODEL_WRITE/SHOT_WRITE/NOWRITE_ONCE \SPECTROSCOPY::TOP.HIREXSR.RESULTS.BSTI:FIT/COMPRESS_ON_PUT/NODO_NOT_COMPRESS
PUT \SPECTROSCOPY::TOP.HIREXSR.RESULTS.BSTI:FIT ""
SET NODE /ON/NOINCLUDE \SPECTROSCOPY::TOP.HIREXSR.RESULTS.BSTI:FIT
SET NODE/MODEL_WRITE/SHOT_WRITE/NOWRITE_ONCE \SPECTROSCOPY::TOP.HIREXSR.RESULTS.BSTI:INST/COMPRESS_ON_PUT/NODO_NOT_COMPRESS
PUT \SPECTROSCOPY::TOP.HIREXSR.RESULTS.BSTI:INST ""
SET NODE /ON/NOINCLUDE \SPECTROSCOPY::TOP.HIREXSR.RESULTS.BSTI:INST
WRITE
