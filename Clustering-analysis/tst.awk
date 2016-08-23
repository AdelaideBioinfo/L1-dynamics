BEGIN { RS=""; ORS="\n\n"; FS="\n"; OFS=", " }
{
    delete keys
    for (i=1; i<=NF; i++) {
        split($i,f," ")
        keys[f[2]]
    }
}
length(keys) > 1 {
    i=0
    keyList=""
    for (key in keys) {
        keyList = keyList (++i>1?OFS:"") key
    }
    if (!seen[keyList]++) {
        print NR, keyList
    }
}
