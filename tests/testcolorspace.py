from cutadapt.colorspace import encode, decode

# If there are any unknown characters in the test sequence,
# round tripping will only work if all characters after the
# first unknown character are also unknown:
# encode("TNGN") == "T444", but
# decode("T444") == "TNNN".

sequences = [
	"",
	"C",
	"ACGGTC",
	"TN",
	"TN.",
	"TNN.N",
	"CCGGCAGCATTCATTACGACAACGTGGCACCGTGTTTTCTCGGTGGTA",
	"TGCAGTTGATGATCGAAGAAAACGACATCATCAGCCAGCAAGTGC",
	"CAGGGTTTGATGAGTGGCTGTGGGTGCTGGCGTATCCGGG"
	]


def test_encode():
	assert encode("AA") == "A0"
	assert encode("AC") == "A1"
	assert encode("AG") == "A2"
	assert encode("AT") == "A3"
	assert encode("CA") == "C1"
	assert encode("CC") == "C0"
	assert encode("CG") == "C3"
	assert encode("CT") == "C2"
	assert encode("GA") == "G2"
	assert encode("GC") == "G3"
	assert encode("GG") == "G0"
	assert encode("GT") == "G1"
	assert encode("TA") == "T3"
	assert encode("TC") == "T2"
	assert encode("TG") == "T1"
	assert encode("TT") == "T0"

	assert encode("TN") == "T4"
	assert encode("NT") == "N4"
	assert encode("NN") == "N4"

	assert encode("ACGGTC") == "A13012"
	assert encode("TTT.N") == "T0044"
	assert encode("TTNT.N") == "T04444"


def test_decode():
	for s in sequences:
		expected = s.replace('.', 'N')
		encoded = encode(s)
		assert decode(encoded) == expected
	assert decode('A.') == 'AN'
	assert decode('C.') == 'CN'
	assert decode('G.') == 'GN'
	assert decode('T.') == 'TN'
