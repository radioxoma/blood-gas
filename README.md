# Acid-base status, blood gas analysis

Some time ago, being a medical student, I started learning acid-base balance in our hospital [ICU](https://en.wikipedia.org/wiki/Intensive_care_unit). We have [one](http://www.radiometer.com/en/products/blood-gas-testing/abl800-flex-blood-gas-analyzer) blood gas analyzer and I always wondering why it report is so complicated. I know only that many of those parameters are derived from few really measured values.

Information about applied calculations relatively scattered; online calculators use simplified techniques. But fortunately there is [ABL800 Flex Reference manual](http://www.radiometeramerica.com/~/media/Files/RadiometerComCloneset/RAME/Manuals/ABL800/989-963I%20ABL800%20Reference%20Manual%20-%20English%20US.pdf), which contains formulas implemented in device software.

Application of mentioned formulas is not always straightforward, so I decided make some reverse-engineering of device calculation function in learning purposes.
