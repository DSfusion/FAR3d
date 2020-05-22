These directories contain typical VMEC and BOOZ_XFORM files. The input.*
files are VMEC inputs and the in_booz.* files are BOOZ_XFORM inputs. The
BOOZ_XFORM inputs typically start with a requested number of toroidal/polodal
modes for the transformed data (XBOOZ_XFORM may increase these as needed
for accuracy). Next, the file extension for the associated VMEC output file
is given. Finally, a list of flux surfaces is provided where the tranform data
is needed. Typically, these include all of the surfaces in the VMEC run, except
the first and last surfaces, since these can sometimes have spurious data.
