#!/usr/bin/env python
import re

def foo():
    d = open("studies_tmp.html").read()

    # Parsing HTML with regular expressions is bad, but because this is
    # generated, this will always work.
    x = re.compile('<dt>.*?"([^"]+?)".*?<dd>(.+?)</dd>', re.DOTALL)

    entries = x.findall(d)
    year = re.compile(".*?-([0-9]{4})")
    res = [[int(year.match(i[0]).groups()[0]), 
            i[1].replace(">DOI<", ">doi<")] for i in entries]

    # This is ugly, but it works (just want the unique years,
    # revsorted).
    years = [i for i in set([i[0] for i in res])]
    years.sort(reverse=True)

    ret = []
    for y in years:
        if y == max(years):
            ret.append('<h2>%d/in press</h2>' % y)
        else:
            ret.append('<h2>%d</h2>' % y)
        ret.append('<ul class="studies">')
        ret += ["<li>\n%s\n</li>" % i[1].strip() for i in res if i[0] == y]
        ret.append("</ul>\n")
        
    str = open("template.html").read() % "\n".join(ret)
    open("studies.html", "w+").write(str)

if __name__ == "__main__":
    foo()
