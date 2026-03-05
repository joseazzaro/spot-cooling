import pdfplumber, re
pdf = 'ASHRAE-D-HO-2668.pdf'
pat = re.compile(r'Example\s*1|Table\s*2|Table\s*3|Table\s*4|Table\s*5', re.I)
with pdfplumber.open(pdf) as d:
    print('PAGES', len(d.pages))
    found = 0
    for i, p in enumerate(d.pages, 1):
        t = p.extract_text() or ''
        if pat.search(t):
            found += 1
            print(f'\n=== PAGE {i} ===')
            for ln in t.splitlines()[:45]:
                print(ln)
    print('\nFOUND', found)
