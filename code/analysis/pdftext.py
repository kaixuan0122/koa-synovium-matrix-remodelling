import re, zlib, pathlib

def _unescape(b):
    out=bytearray(); i=0
    while i < len(b):
        c=b[i]
        if c==0x5c and i+1 < len(b):          # backslash
            n=b[i+1]
            if n in b'nrtbf': out.append({0x6e:10,0x72:13,0x74:9,0x62:8,0x66:12}[n])
            elif 0x30 <= n <= 0x37:
                j=i+1; oct_=b''
                while j < len(b) and len(oct_)<3 and 0x30 <= b[j] <= 0x37: oct_+=bytes([b[j]]); j+=1
                out.append(int(oct_,8)); i=j; continue
            else: out.append(n)
            i+=2; continue
        out.append(c); i+=1
    return bytes(out).decode('latin-1')

def strings_in(seg):
    """按 PDF 语法取字符串字面量，正确处理 \\( \\) \\\\ 等转义。"""
    res=[]; i=0
    while i < len(seg):
        if seg[i:i+1]==b'(':
            depth=1; j=i+1; buf=bytearray()
            while j < len(seg) and depth>0:
                c=seg[j]
                if c==0x5c: buf+=seg[j:j+2]; j+=2; continue
                if c==0x28: depth+=1
                elif c==0x29:
                    depth-=1
                    if depth==0: break
                buf.append(c); j+=1
            res.append(_unescape(bytes(buf))); i=j+1; continue
        i+=1
    return res

def pdf_text(path):
    d=pathlib.Path(path).read_bytes(); recs=[]
    for m in re.finditer(rb'stream\r?\n', d):
        s=m.end(); e=d.find(b'endstream', s)
        if e<0: continue
        try: raw=zlib.decompress(d[s:e])
        except Exception: continue
        for bm in re.finditer(rb'BT(.*?)ET', raw, re.S):
            seg=bm.group(1); x=y=0.0
            for tm in re.finditer(rb'([-\d.]+) ([-\d.]+) ([-\d.]+) ([-\d.]+) ([-\d.]+) ([-\d.]+) Tm', seg):
                x=float(tm.group(5)); y=float(tm.group(6))
            for item in re.finditer(rb'(\[[^\]]*\]|\([^)]*\)|\((?:\\.|[^)])*\))\s*(Tj|TJ)', seg):
                t=''.join(strings_in(item.group(1)))
                if t.strip(): recs.append((round(y,1), round(x,1), t))
    return recs, d
