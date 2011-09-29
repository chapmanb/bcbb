"""Code for writing to google docs
"""
def _from_unicode(unistr,encoding='utf-8'):
    """Encode a unicode string, by default into utf-8"""
    if isinstance(unistr,unicode):
        unistr = unistr.encode(encoding)
    return unistr

def _to_unicode(str,encoding='utf-8'):
    """Decode a string into unicode"""
    if isinstance(str,basestring):
        if not isinstance(str,unicode):
            str = unicode(str,encoding)
    return str
    