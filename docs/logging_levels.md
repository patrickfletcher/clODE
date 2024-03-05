# Logging levels

clODE has several internal logging levels. These are set by the
`clode.set_log_level` function, which takes an enum argument. The
levels are:

* `clode.log_level.off`: No logging
* `clode.log_level.critical`: Only log critical errors
* `clode.log_level.error`: Only log errors
* `clode.log_level.warning`: Log errors and warnings
* `clode.log_level.info`: Log errors, warnings, and info
* `clode.log_level.debug`: Log errors, warnings, info, and debug
* `clode.log_level.trace`: Log everything, including trace messages

To set the logging level, use the `clode.set_log_level` function:

```python
import clode

clode.set_log_level(clode.log_level.debug) # Set the logging level to debug

# Do some stuff

clode.set_log_level(clode.log_level.off) # Turn off logging
```
