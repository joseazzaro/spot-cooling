import pdfplumber
pdf='ASHRAE-D-HO-2668.pdf'
with pdfplumber.open(pdf) as d:
    for i in [5,6,7,8,9,10,11,12]:
        if i>len(d.pages):
            continue
        p=d.pages[i-1]
        tabs=p.extract_tables()
        if not tabs:
            continue
        print(f'\n=== PAGE {i} TABLES {len(tabs)} ===')
        for ti,t in enumerate(tabs,1):
            print(f'-- Table {ti} rows={len(t)} cols={max((len(r) for r in t if r), default=0)}')
            for r in t[:12]:
                print(r)
