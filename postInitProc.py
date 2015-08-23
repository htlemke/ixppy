def process(d):
  if hasattr(d,'evrBool'):
    try:
      d['eventCodeBool'] = [te==1 for te in d.evrBool.data.get_memdata()]
    except:
      print "Post process unpacking evrBool data did not succeed for some reason!!!"

  if hasattr(d,'adc'):
    try:
      d['adcV'] = d.evrBool.data.get_memdata()


    except:
      print "Post process unpacking evrBool data did not succeed for some reason!!!"

