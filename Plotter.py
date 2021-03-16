#Created 12.03.2020 by Ole Oggesen

import numpy as np
import matplotlib.pyplot as plt
import os
import array
import xml.etree.ElementTree as ET


def plot(xarray, yarray, title, xlabel, ylabel,  save, labelarray = None, xlim = None, ylim = None):
    '''Function for plotting the measurments on the NMR MRT
    y can contain multiple dataarrays,
    save is either True or False
    xlim is a touple with 2 values
    Needs a plots folder in same folder to save the plots'''
    plt.figure()
    for x, y in zip(xarray, yarray):
        plt.plot(x, y, linewidth=1)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if labelarray:
        plt.legend(labelarray)
    if xlim is not None:
        if isinstance(xlim, type(())):
            plt.xlim(xlim)
    if ylim is not None:
        if isinstance(ylim, type(())):
            plt.ylim(ylim)
    if save == True:
        savetitle = "./plots/" + title.replace(" ", "_") + ".pdf"
        plt.savefig(savetitle, format="pdf")
        plt.close()
    else:
        plt.show()

def load(path):
    '''Loads data
    supports datatypes .bin and .csv'''
    if ".bin" in path:
        with open(path, "rb") as file:
            data = file.read()
            if data[0:5] == b'<?xml':
                tree = ET.parse(path)
                data = tree.getroot()
            else:
                data = array.array("f", data)
            return data
    elif ".csv" in path:
        with open(path, 'r') as f:
            data = f.read()
            data = data.split("\n")
            try:
                data.remove("")
            finally:
                assert isinstance(data, list)
                return data

def parseData(path):
    """loads the data from the .Wfm file
     and parses all data needed for plotting
      in args from the casual file"""

    try:
        data = load(path + ".Wfm.csv")
        args = load(path + ".csv")
    except FileNotFoundError:
        try:
            data = load(path + ".Wfm.bin")
            args = load(path + ".bin")
        except FileNotFoundError:
            raise FileNotFoundError("File not found. Only bin and csv datatypes are supported.")
    val = {}

    # parse for x values
    if isinstance(args, type([])):
        for arg in args:
            if "HardwareXStart" in arg:
                tmp = arg.split(":")
                val["HardwareXStart"] = float(tmp[1])
            elif "HardwareXStop" in arg:
                tmp = arg.split(":")
                val["HardwareXStop"] = float(tmp[1])

        #parse for y values
        for i, point in enumerate(data):
            data[i] = float(point)
        data = [data]  # TODO check if .csv files with more then one waveform are read correctly

    if isinstance(args, ET.Element):
        for child in args:
            for kids in child:
                attribute = kids.attrib
                if "HardwareXStart" in attribute["Name"]:
                    val["HardwareXStart"] = float(attribute["Value"])
                if "HardwareXStop" in attribute["Name"]:
                    val["HardwareXStop"] = float(attribute["Value"])
                if "SignalHardwareRecordLength" in attribute["Name"]:
                    #decodes how many signals were recorded and seperates them
                    numOfSignals = int(len(data)/int(attribute["Value"]))
                    tmp = []
                    for i in range(numOfSignals):
                        tmp.append(data[i::numOfSignals])
                    data = tmp

    assert len(val) == 2
    return data,  val


def Plotter(paths, title, save=False, labelarray = None, xlim = None, extraData = None):
    """Plots all data from given path automaticly
    extraData for extra plots will always be plotted at last"""
    yarray = []
    xarray = []

    for path in paths:
        data, args = parseData(path)

        for dataElem in data:
            if xlim is not None:
                dataElem, x = xWindow(dataElem, xlim, args)
            else:
                x = np.linspace(args["HardwareXStart"], args["HardwareXStop"], len(dataElem))
            if save is True:
                maxlength = 3000
            else:
                maxlength = 10000
            if len(dataElem) > maxlength:
                x, dataElem = downSample(data=dataElem, maxlength=maxlength, x=x)
            yarray.append(list(dataElem))
            xarray.append(list(x))
    if extraData is not None:
        xarray.append(list(extraData[0]))
        yarray.append(list(extraData[1]))
    counter = 0
    for i in range(len(xarray)):
        while np.max(xarray[i]) < 0.1 and np.min(xarray[i]) > -0.1:
            for k in range(len(xarray)):
                for j in range(len(xarray[k])):
                    xarray[k][j] = xarray[k][j] * 1000
            if xlim is not None:
                xlim = (xlim[0] * 1000, xlim[1] * 1000)
            counter = counter + 1
        xlabellist = ["s", "ms", "us", "ns", "ps"]
        xlabel = "time in " + xlabellist[counter]

    try:
        assert len(xarray) != 0 and len(yarray) != 0
    except:
        raise AssertionError("This file in contains no information. An empty list was created.")

    plot(xarray, yarray, title, xlabel, "voltage in V", save, labelarray=labelarray, xlim=xlim)

def findData(path):
    paths = []
    with os.scandir(path) as it:
        for entry in it:
            if entry.is_dir():
                try:
                    paths.extend(findData(entry.path))
                except :
                    pass
            elif entry.is_file():
                if ".Wfm.csv" in entry.name:
                    tmp = str(entry.path).replace(".Wfm.csv", "")
                    paths.append(tmp)
                elif ".Wfm.bin" in entry.name:
                    tmp = str(entry.path).replace(".Wfm.bin", "")
                    paths.append(tmp)
    if paths == []:
        pass
    else:
        return paths


def Autoplot(path, titledict = None):
    paths = findData(path)
    for i, path in enumerate(paths):
        title = None
        if titledict is not None:
            for key in titledict.keys():
                if key in path:
                    title = titledict[key]
                    break
        if title == None:
            title = path.split("/")
            title = title[-1]
        Plotter([path], title, True)


def FFT(path, title=None, xlim=None, save=False, freqLim=None):
    data, args = parseData(path)

    if xlim is not None:
        data, x = xWindow(data[0], xlim, args)
        x = None
        step = (xlim[1]-xlim[0])/len(data)
    else:
        data = data[0]
        step = (args["HardwareXStop"] - args["HardwareXStart"]) / len(data)

    Tmax=0.25e-7    #1/fs
    if step < Tmax:
        maxlength = int(step/Tmax*len(data))
        data = downSample(data=data, maxlength=maxlength)
        step = int(Tmax/step)*step

    N = len(data)
    fft = np.fft.fft(data, N)
    data = None
    mag = 2*abs(fft)/N
    fft = None
    freq = np.fft.fftfreq(N, step)
    for i in range(len(freq)):
        freq[i] = freq[i]*1e-6
    freqLim = (freqLim[0]*1e-6, freqLim[1]*1e-6)


    plt.figure()
    plt.stem(freq, mag, use_line_collection=True)
    plt.xlabel("frequency in MHz")
    plt.ylabel("amplitude in V")
    if freqLim is not None:
        plt.xlim(freqLim)
    if title is not None:
        plt.title(title)
    if save is True:
        if title is not None:
            saveTitle = "./plots/" + title.replace(" ", "_") + "_FFT.pdf"
            plt.savefig(saveTitle, format="pdf")
        else:
            raise NotImplementedError("Title is needed to save the figure!")
        plt.close()
    else:
        plt.show()


def xWindow(dataElem, xlim, args):
    step = (args["HardwareXStop"] - args["HardwareXStart"]) / len(dataElem)
    dataElem = dataElem[int((xlim[0] - args["HardwareXStart"]) / step):int((xlim[1] - args["HardwareXStart"]) / step)]
    x = np.linspace(xlim[0], xlim[1], len(dataElem))
    try:
        assert len(dataElem) != 0
    except AssertionError:
        raise AssertionError("Wrong xlim data input. There is no data for given time period.")
    return dataElem, x

def downSample(data, maxlength, x=None):
    data = data[0::int(len(data) / maxlength)]
    if x is not None:
        x = x[0::int(len(x) / maxlength)]
        return x, data
    return data

def getData(path, xlim=None):
    data, args = parseData(path)

    if xlim is not None:
        data, x = xWindow(data[0], xlim, args)
    return data

def calcPower(data):
    #Calculates power over a 50 Ohms resistor
    return np.sum(np.array(data)**2)/50/len(data)










if __name__ == '__main__':
    #TODO titledict only works if titles are keys are unique for every path
    """
    titledict = {}
    titledict["GateAndRFPulseAmplified0212"] = "RF and gate pulse amplified"
    titledict["GateAndRFPulseNotAmplified0212"] = "RF and gate pulse not amplified"
    titledict["90und180degPulse50Ohm"] = "90 and 180 deg  RF pulse"
    titledict["dummyCoilGateRinging"] = "50 Ohms Dummy Coil Ringing"
    titledict["NoCoilNoTransSignalGateRinging"] = "Dummy Coil No Tx Signal Ringing"
    titledict["RealCoilGateRinging"] = "Attached RF Coil Ringing"
    titledict["TRSwitchJustGatePulse"] = "Gate Ringing with no TX Signal"
    titledict["RFPulsePowerAmplifier90and180deg"] = "Output RF Amplifier"
    titledict["90degPulse"] = "90 degree pulse"

    curr_path = os.path.dirname(__file__)   #returns the current path of the python skript
    curr_path = os.path.dirname(curr_path)  # ".."

    Autoplot(curr_path, titledict = titledict)
    """

    save = True

    """
    start = 0.000025
    end = 0.000026
    xlim = (start, end)
    t = np.linspace(start, end, 3000)
    y = 0.15*np.sin(2*np.pi*4.3576e6*t+0.39*np.pi)
    xData = [t, y]
    path = "/home/ole/StorageServer/Work/BachelorThesis/RedPitaya16BitUndMessaufbau/Messergebnisse/11.01.2021/90degPulse"
    Plotter(paths=[path], title="Output RP", save=save, xlim=xlim, labelarray=["output", "ideal"], extraData=xData)

    paths = ["/home/ole/StorageServer/Work/BachelorThesis/RedPitaya16BitUndMessaufbau/Messergebnisse/11.01.2021/90and180degPulse", "/home/ole/StorageServer/Work/BachelorThesis/RedPitaya16BitUndMessaufbau/Messergebnisse/11.01.2021/90and180degPulse50Ohm"]

    Plotter(paths=paths, title="Difference in input resistance", save=save, labelarray=["1 MOhm", "50 Ohm"])


    freqLim = (3e6, 6e6)
    path = "/home/ole/StorageServer/Work/BachelorThesis/RedPitaya16BitUndMessaufbau/Messergebnisse/11.01.2021/90degPulse"
    xlim = (0.0, 0.00003)

    FFT(path=path, title="Output RP", xlim=xlim, save=save, freqLim=freqLim)

    
    path = "/home/ole/StorageServer/Work/BachelorThesis/RedPitaya16BitUndMessaufbau/Messergebnisse/13.01.2021/TRSwitchJustGatePulse"

    xlim = (-0.0001, 0.0001)
    freqLim = (0, 5e6)

    Plotter(paths=[path], title="Gate Ringing Zoom", save=save, xlim=xlim)

    xlim = (0, 0.0001)
    FFT(path=path, title="Gate Ringing TR Switch", xlim=xlim, freqLim=freqLim, save=save)

    #calculate power of pulse with 50 Ohms
    voltage = getData(path=path, xlim=(0, 40e-6))
    power = calcPower(data=voltage)
    print("Power of the first pulse: {} W".format(power))

    freqLim = (0.0, 5e6)
    path = "/home/ole/StorageServer/Work/BachelorThesis/RedPitaya16BitUndMessaufbau/Messergebnisse/13.01.2021/RFPulsePowerAmplifier90and180deg"
    xlim=(-0.00015, 0.0001)
    Plotter(paths=[path], title="90 deg pulse power amplifier", save=save, xlim=xlim)
    """
    input = np.array([i*0.1 for i in range(1, 13)])
    theo_output = np.array(input)*1000
    real_output = np.array([98.98, 191.3192, 281.7609, 373.1514, 461.6958, 559.7271, 632.46, 711.5175, 787.4127, 860.1456, 936.0408, 1011.936])

    plot([input, input], [theo_output, real_output], "RF amplifier", "voltage in V", "voltage in V", True, ["ideal", "measured"], (0.1, 0.6), (100, 600))


