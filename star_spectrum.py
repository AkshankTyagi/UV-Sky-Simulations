def GET_STAR_TEMP(sptype):
    sptype = hipline["sp_type"]

    if sptype.startswith("sd"):
        sptype = sptype[2:]

    if sptype.startswith("O3"):
        hipline["temperature"] = 5
    elif sptype.startswith("O4"):
        hipline["temperature"] = 7
    elif sptype.startswith("O5"):
        hipline["temperature"] = 10
    elif sptype.startswith("O6"):
        hipline["temperature"] = 11
    elif sptype.startswith("O7"):
        hipline["temperature"] = 13
    elif sptype.startswith("O8"):
        hipline["temperature"] = 15
    elif sptype.startswith("O9"):
        hipline["temperature"] = 18
    elif sptype.startswith("B0"):
        hipline["temperature"] = 20
    elif sptype.startswith("B1"):
        hipline["temperature"] = 25
    elif sptype.startswith("B2"):
        hipline["temperature"] = 28
    elif sptype.startswith("B3"):
        hipline["temperature"] = 31
    elif sptype.startswith("B4"):
        hipline["temperature"] = 33
    elif sptype.startswith("B5"):
        hipline["temperature"] = 35
    elif sptype.startswith("B6"):
        hipline["temperature"] = 36
    elif sptype.startswith("B7"):
        hipline["temperature"] = 37
    elif sptype.startswith("B8"):
        hipline["temperature"] = 41
    elif sptype.startswith("B9"):
        hipline["temperature"] = 49
    elif sptype.startswith("A0"):
        hipline["temperature"] = 51
    elif sptype.startswith("A1"):
        hipline["temperature"] = 52
    elif sptype.startswith("A2"):
        hipline["temperature"] = 53
    elif sptype.startswith("A3"):
        hipline["temperature"] = 56
    elif sptype.startswith("A4"):
        hipline["temperature"] = 56
    elif sptype.startswith("A5"):
        hipline["temperature"] = 56
    elif sptype.startswith("A6"):
        hipline["temperature"] = 57
    elif sptype.startswith("A7"):
        hipline["temperature"] = 58
    elif sptype.startswith("A8"):
        hipline["temperature"] = 59
    elif sptype.startswith("A9"):
        hipline["temperature"] = 60
    elif sptype.startswith("F0"):
        hipline["temperature"] = 60
    elif sptype.startswith("F1"):
        hipline["temperature"] = 63
    elif sptype.startswith("F2"):
        hipline["temperature"] = 61
    elif sptype.startswith("F3"):
        hipline["temperature"] = 62
    elif sptype.startswith("F4"):
        hipline["temperature"] = 62
    elif sptype.startswith("F5"):
        hipline["temperature"] = 63
    elif sptype.startswith("F6"):
        hipline["temperature"] = 63
    elif sptype.startswith("F7"):
        hipline["temperature"] = 63
    elif sptype.startswith("F8"):
        hipline["temperature"] = 64
    elif sptype.startswith("F9"):
        hipline["temperature"] = 65
    elif sptype.startswith("G0"):
        hipline["temperature"] = 65
    elif sptype.startswith("G1"):
        hipline["temperature"] = 66
    elif sptype.startswith("G2"):
        hipline["temperature"] = 66
    elif sptype.startswith("G3"):
        hipline["temperature"] = 66
    elif sptype.startswith("G4"):
        hipline["temperature"] = 66
    elif sptype.startswith("G5"):
        hipline["temperature"] = 66
    elif sptype.startswith("G6"):
        hipline["temperature"] = 66
    elif sptype.startswith("G7"):
        hipline["temperature"] = 67
    elif sptype.startswith("G8"):
        hipline["temperature"] = 67
    elif sptype.startswith("G9"):
        hipline["temperature"] = 67
    elif sptype.startswith("K0"):
        hipline["temperature"] = 68
    elif sptype.startswith("K1"):
        hipline["temperature"] = 69
    elif sptype.startswith("K2"):
        hipline["temperature"] = 70
    elif sptype.startswith("K3"):
        hipline["temperature"] = 70
    elif sptype.startswith("K4"):
        hipline["temperature"] = 71
    elif sptype.startswith("K5"):
        hipline["temperature"] = 72
    elif sptype.startswith("K6"):
        hipline["temperature"] = 72
    elif sptype.startswith("K7"):
        hipline["temperature"] = 73
    elif sptype.startswith("K8"):
        hipline["temperature"] = 73
    elif sptype.startswith("K9"):
        hipline["temperature"] = 73
    elif sptype.startswith("M0"):
        hipline["temperature"] = 74
    elif sptype.startswith("M1"):
        hipline["temperature"] = 74
    elif sptype.startswith("M2"):
        hipline["temperature"] = 75
    elif sptype.startswith("M3"):
        hipline["temperature"] = 75
    elif sptype.startswith("M4"):
        hipline["temperature"] = 75
    elif sptype.startswith("M5"):
        hipline["temperature"] = 75
    elif sptype.startswith("M6"):
        hipline["temperature"] = 75
    elif sptype.startswith("M7"):
        hipline["temperature"] = 75
    elif sptype.startswith("M8"):
        hipline["temperature"] = 75
    elif sptype.startswith("M9"):
        hipline["temperature"] = 75
    elif sptype.startswith("O"):
        hipline["temperature"] = 13
    elif sptype.startswith("B"):
        hipline["temperature"] = 35
    elif sptype.startswith("M"):
        hipline["temperature"] = 75
    elif sptype.startswith("C"):
        hipline["temperature"] = 75
    elif sptype.startswith("A"):
        hipline["temperature"] = 56
    elif sptype.startswith("R"):
        hipline["temperature"] = 75
    elif sptype.startswith("G"):
        hipline["temperature"] = 66
    elif sptype.startswith("W"):
        hipline["temperature"] = 35
    elif sptype.startswith("K"):
        hipline["temperature"] = 72
    elif sptype.startswith("N"):
        hipline["temperature"] = 72
    elif sptype.startswith("S"):
        hipline["temperature"] = 72
    elif sptype.startswith("F"):
        hipline["temperature"] = 63
    elif sptype.startswith("DA"):
        hipline["temperature"] = 35
    else:
        # Handle other spectral types or unknown cases
        hipline["temperature"] = 66  # Default temperature

    return hipline["temperature"]

# Example usage
hipline = {"sp_type": "sdB3V", "temperature": 0}
temperature = GET_STAR_TEMP(hipline)
print(f"Temperature: {temperature}")
