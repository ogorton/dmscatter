f = open("validation.md","w+")
f.write("---\ngeometry:\n- margin=.5in\n---\n# DMFortFactor Validation Plots for $^{29}$Si\n")
pairs = [
        (1,1), (1,2), (1,3),
        (2,2), (2,3),
        (3,3),
        (4,4), (4,5),
        (5,5), (5,6),
        (6,6),
        (7,7),
        (8,8), (8,9),
        (9,9),
        (10,10),
        (11,11), (11,12), (11,15),
        (12,12), (12,15),
        (13,13),
        (14,14),
        (15,15)]

for pair in pairs:
    o1, o2 = pair[:]
    f.write("![](log.usd.n.%i.%i.pdf){width=50%%}![](perr.log.usd.n.%i.%i.pdf){width=50%%}\n"%(o1,o2,o1,o2))
    f.write("![](log.usd.p.%i.%i.pdf){width=50%%}![](perr.log.usd.p.%i.%i.pdf){width=50%%}\n"%(o1,o2,o1,o2))

f.close()    
    
