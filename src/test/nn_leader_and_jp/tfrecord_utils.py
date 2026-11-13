"""Minimal dependency-free TFRecord/protobuf fixtures for CLI tests."""
import struct


def masked_crc(data):
    # Standard reflected CRC32C and TFRecord mask, without extra dependencies.
    crc = 0xffffffff
    for byte in data:
        crc ^= byte
        for _ in range(8):
            crc = (crc >> 1) ^ (0x82f63b78 if crc & 1 else 0)
    crc ^= 0xffffffff
    return struct.pack('<I', (((crc >> 15) | (crc << 17)) + 0xa282ead8) & 0xffffffff)


def string_field(number, value):
    data = value.encode()
    assert len(data) < 128
    return bytes([number * 8 + 2, len(data)]) + data


def tfrecord(name, neighbours=()):
    # nnbr.NearNeighbours: smiles = 1, name = 2, repeated Nbr = 5.
    proto = string_field(1, 'C') + string_field(2, name)
    for neighbour, distance in neighbours:
        nbr = string_field(2, neighbour) + b'\x1d' + struct.pack('<f', distance)
        proto += bytes([42, len(nbr)]) + nbr
    return frame(proto)



def frame(proto):
    length = struct.pack('<Q', len(proto))
    return length + masked_crc(length) + proto + masked_crc(proto)
