# Logging levels

clODE exposes logging levels through the `LogLevel` enum and the `set_log_level(...)` helper.

Available levels:

- `clode.LogLevel.off`: no logging
- `clode.LogLevel.critical`: only critical failures
- `clode.LogLevel.err`: errors and above
- `clode.LogLevel.warn`: warnings and above
- `clode.LogLevel.info`: informational messages and above
- `clode.LogLevel.debug`: debug output and above
- `clode.LogLevel.trace`: the most verbose logging

To change the active level:

```python
import clode

clode.set_log_level(clode.LogLevel.debug)

# Do some stuff

clode.set_log_level(clode.LogLevel.off)
```

You can inspect the current level with:

```python
current_level = clode.get_log_level()
print(current_level)
```
