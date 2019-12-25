# Acid-base status, blood gas analysis

> This code partially used in my clinical application - [Heval](https://github.com/radioxoma/heval). Feel free to check it out.

Some time ago, being a medical student, I started learning acid-base balance in our hospital [ICU](https://en.wikipedia.org/wiki/Intensive_care_unit). We been using [ABL 800 Flex](https://www.radiometer.com/en/products/blood-gas-testing/abl800-flex-blood-gas-analyzer) blood gas analyzer and I always wondering why it's report is so complicated. I know only that many of those parameters are derived from few really measured values.

Information about applied calculations relatively scattered; online calculators use simplified techniques. But fortunately there is [ABL800 Flex Reference manual](http://www.radiometeramerica.com/~/media/Files/RadiometerComCloneset/RAME/Manuals/ABL800/989-963I%20ABL800%20Reference%20Manual%20-%20English%20US.pdf), which contains formulas implemented in device software.

Application of mentioned formulas is not always straightforward, so I decided to make some reverse-engineering of device calculation function in edication purposes.
