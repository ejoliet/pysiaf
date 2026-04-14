from astropy.table import Table

oss_data = {
    "InstrName": [],
    "AperName": [],
    "DDCName": [],
    "OSS_Version": [],
    "Action": [],
    "V2Ref": [],
    "V3Ref": [],
    "V3IdlYAngle": []
}

all_data = []
data = [
    "FGS",
    "*",
    "*",
    "8.4",
    "default",
    "default",
    "default",
    "default"
]
all_data.append(data)
data = [
    "FGS",
    "J-FRAME",
    "*",
    "8.4",
    "default",
    180.255,
    175.783,
    0.04168887
]
all_data.append(data)
data = [
    "FGS",
    "J-FRAME",
    "*",
    "10.1",
    "append",
    196.3,
    181.806,
    0.05134242
]
all_data.append(data)
data = [
    "FGS",
    "J-FRAME",
    "*",
    "11.1",
    "append",
    205.895,
    76.798,
    0.0513463,
]
all_data.append(data)
data = [
    "MIRI",
    "*",
    "*",
    "8.4",
    "default",
    "default",
    "default",
    "default"
]
all_data.append(data)
data = [
    "NIRCAM",
    "*",
    "*",
    "8.4",
    "default",
    "default",
    "default",
    "default"
]
all_data.append(data)
data = [
    "NIRCAM",
    "NRCA3_SUB64P",
    "*",
    "11.1",
    "default",
    "default",
    "default",
    "default"
]
all_data.append(data)
data = [
    "NIRCAM",
    "NRCA4_SUB160P",
    "*",
    "11.1",
    "default",
    "default",
    "default",
    "default"
]
all_data.append(data)
data = [
    "NIRCAM",
    "NRCA4_SUB400P",
    "*",
    "11.1",
    "default",
    "default",
    "default",
    "default"
]
all_data.append(data)
data = [
    "NIRCAM",
    "NRCA5_SUB64P",
    "NRCA_CNTR",
    "8.4",
    "default",
    "default",
    "default",
    "default"
]
all_data.append(data)
data = [
    "NIRCAM",
    "NRCA5_SUB64P",
    "NRC_MASKLWB",
    "11.1",
    "default",
    "default",
    "default",
    "default"
]
all_data.append(data)
data = [
    "NIRCAM",
    "NRCA5_SUB160P",
    "NRCALL_CNTR",
    "8.4",
    "default",
    "default",
    "default",
    "default"
]
all_data.append(data)
data = [
    "NIRCAM",
    "NRCA5_SUB160P",
    "NRC_MASKLWB",
    "11.1",
    "default",
    "default",
    "default",
    "default"
]
all_data.append(data)
data = [
    "NIRCAM",
    "NRCA5_SUB400P",
    "NRCALL_CNTR",
    "8.4",
    "default",
    "default",
    "default",
    "default"
]
all_data.append(data)
data = [
    "NIRCAM",
    "NRCA5_SUB400P",
    "NRC_MASKLWB",
    "11.1",
    "default",
    "default",
    "default",
    "default"
]
all_data.append(data)
data = [
    "NIRISS",
    "*",
    "*",
    "8.4",
    "default",
    "default",
    "default",
    "default"
]
all_data.append(data)
data = [
    "NIRSPEC",
    "*",
    "*",
    "8.4",
    "default",
    "default",
    "default",
    "default"
]
all_data.append(data)


for row in all_data:
    oss_data["InstrName"].append(row[0])
    oss_data["AperName"].append(row[1])
    oss_data["DDCName"].append(row[2])
    oss_data["OSS_Version"].append(row[3])
    oss_data["Action"].append(row[4])
    oss_data["V2Ref"].append(row[5])
    oss_data["V3Ref"].append(row[6])
    oss_data["V3IdlYAngle"].append(row[7])

t = Table(oss_data)
t.write("OSS_VERSION_TABLE.txt", format="ascii.fixed_width", overwrite=True)
