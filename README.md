# GTA VI Twitch CCU forecast

[![make.fun · gta6-twitch-ccu-forecast-01](https://make.fun/badge/gta6-twitch-ccu-forecast-01.svg)](https://make.fun/bounty/gta6-twitch-ccu-forecast-01?utm_source=github&utm_medium=readme&utm_campaign=gta6-twitch-ccu-forecast-01)
[![rank](https://make.fun/badge/gta6-twitch-ccu-forecast-01/Amburadullzz/gta6-ccu-forecast/rank.svg)](https://make.fun/bounty/gta6-twitch-ccu-forecast-01/leaderboard?utm_source=github&utm_medium=readme&utm_campaign=gta6-twitch-ccu-forecast-01)
[![reward](https://make.fun/badge/gta6-twitch-ccu-forecast-01/Amburadullzz/gta6-ccu-forecast/reward.svg)](https://make.fun/bounty/gta6-twitch-ccu-forecast-01?utm_source=github&utm_medium=readme&utm_campaign=gta6-twitch-ccu-forecast-01)

Hour by hour, for the 24 UTC hours of Grand Theft Auto VI's release day (2026-11-19), five quantiles
of the concurrent Twitch viewers in the game's category, from the day's first hour and the category's
own history.

## The model

Hour 1's print is given, and so is a library of comparable launches: past console releases with their
own hour-by-hour Twitch curves. For every comparable the program measures

    r_c(h) = ln(v_c(h) / v_c(1))

the log ratio of hour `h` to that launch's own first hour. The ratio is scale free, so the shape of a
launch carries across launches of very different size. The ratios are pooled per hour with a Gaussian
kernel on `log10` of hour-1 viewers centred on the case's own first hour, so comparables of a similar
size dominate, and the weighted median is taken so one blow-up curve cannot drag the forecast.

A pooled weighted least-squares slope of the ratio on `log10` of hour-1 viewers, fitted over hours 2
to 8, tilts the curve for a case bigger or smaller than the comparable set. Its input is clipped to the
set's own range, so an unprecedented launch never extrapolates the tilt, and a case beyond every
comparable gets a wider band on top, where no comparable can speak for it. The spread widens with
`sqrt(h - 1)`.

Tuning was leave-one-out over the kit's 37 comparable launches: on the sixteen launches above 100,000
viewers in hour 1 the mean pinball loss is 0.216 against 0.261 for the house baseline's persistence
with a widening band; over all 37 it is 0.421 against 0.930.

## Run it

```
MAKEFUN_INPUT_DIR=input MAKEFUN_OUTPUT_DIR=out python3 main.py
```

Python 3 standard library only, no network, deterministic: two runs over the same case write the same
bytes. One case takes under a second.

Bounty: https://make.fun/bounty/gta6-twitch-ccu-forecast-01
