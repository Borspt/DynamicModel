from Elements import BlockValve

object = BlockValve(elementId=0, status=1, kvFile='regulatorKvFile.csv')

object.calcFlow(openPercentage=75)