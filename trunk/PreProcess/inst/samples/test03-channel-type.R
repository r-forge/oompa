# channel-type-test.ssc

x <- ChannelType('Affymetrix', 'oligo', 100, 100, 'fluor')
print(x)

summary(x)

y <- setDesign(x, 'fake.design')
print(y)
summary(y)
d <- getDesign(y)
d

rm(d, x, y)
