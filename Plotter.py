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
      in args from the casual file. Only compatible with .csv and .bin files"""

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
    """Plots all data from given path automatically, which are stored in .csv or .bin format.
    extraData for extra plots will always be plotted at last"""
    yarray = []
    xarray = []

    for path in paths:
        data, args = parseData(path)

        for dataElem in data:
            if xlim is not None:
                dataElem, x, _ = RectWindow(dataElem, xlim, args)
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

def findData(path, search_key = None):
    """Finds the path of data in a given folder and its subfolders. If search_key is given, only paths which contain the search_key are returned.
    Else it returns all Wfm.bin and Wfm.csv files and excludes its endings."""
    paths = []
    with os.scandir(path) as it:
        for entry in it:
            if entry.is_dir():
                try:
                    paths.extend(findData(entry.path))
                except :
                    pass
            elif entry.is_file():
                if search_key is not None:
                    if search_key in entry.name:
                        paths.append(str(entry.path))
                else:
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
    """Plots all .bin and .csv saved data in given folder and all subfolders."""
    paths = findData(path)
    for path in paths:
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


def FFTPlot(path, title=None, xlim=None, save=False, freqLim=None):
    """Plots a FFT from the data in given path."""
    data, args = parseData(path)

    if xlim is not None:
        data, x, step = RectWindow(data[0], xlim, args)
        x = None
    else:
        data = data[0]  #TODO adapt for multiinput signals
        step = (args["HardwareXStop"] - args["HardwareXStart"]) / len(data)
    Tmax=1e-7

    mag, freq = FFT(data, step, Tmax)
    data = None

    for i in range(len(freq)):
        freq[i] = freq[i] * 1e-6
    freqLim = (freqLim[0] * 1e-6, freqLim[1] * 1e-6)
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


def RectWindow(dataElem, xlim, args):
    step = (args["HardwareXStop"] - args["HardwareXStart"]) / len(dataElem)
    dataElem = dataElem[int((xlim[0] - args["HardwareXStart"]) / step):int((xlim[1] - args["HardwareXStart"]) / step)]
    x = np.linspace(xlim[0], xlim[1], len(dataElem))
    try:
        assert len(dataElem) != 0
    except AssertionError:
        raise AssertionError("Wrong xlim data input. There is no data for given time period.")
    step = (xlim[1] - xlim[0]) / len(dataElem)
    return dataElem, x, step

def downSample(data, maxlength, x=None):
    data = data[0::int(len(data) / maxlength)]
    if x is not None:
        x = x[0::int(len(x) / maxlength)]
        return x, data
    return data

def getData(path, xlim=None):
    data, args = parseData(path)

    if xlim is not None:
        data, x, _ = RectWindow(data[0], xlim, args)
    return data

def calcPower(data):
    #Calculates power over a 50 Ohms resistor
    return np.sum(np.array(data)**2)/50/len(data)

def measureFreq(path, xlim=None):
    """Makes a FFT of the received data, and writes the dominant frequency to an .txt file."""
    data, args = parseData(path)

    for i, d in enumerate(data):
        if xlim is not None:
            d, t, step = RectWindow(dataElem=d, xlim=xlim, args=args)
        else:
            step = (args["HardwareXStop"] - args["HardwareXStart"]) / len(d)
            T0 = step* len(d)
            T0min = 1/10
        if 1/T0 > T0min:
            d.extend(np.zeros(int(T0min / step) - int(T0 / step)))      #Zero padding
        mag, freq = FFT(data=d, step=step)
        freq = freq[0:int(len(freq)/2)]
        idx_start = binarySearch(freq, 3e6)
        idx_stop = binarySearch(freq, 5e6)

        m_max = 0.0
        for f, m in zip(freq[idx_start:idx_stop], mag[idx_start:idx_stop]):
            if m>m_max:
                f_max = f
                m_max = m

        path = os.path.dirname(path)
        with open(path + "/freqChannel{}.txt".format(i+1), "a") as file:
            file.write("{}\n".format(f_max))

def measureFreqAuto(path, xlim = None):
    """Measures the dominant frequency in data stored in path address and saves it to a .txt file."""
    paths = findData(path)
    for path in paths:  #Removes the stored data in .txt files if they exist
        path = os.path.dirname(path)
        for i in range(2):
            data = path + "/freqChannel{}.txt".format(i + 1)
            if os.path.exists(data):
                os.remove(data)
    for path in paths:
        if xlim is not None:
            measureFreq(path=path, xlim = xlim)
        else:
            measureFreq(path=path)
    new_paths = []
    for path in paths:
        path = os.path.dirname(path)
        if path in new_paths:
            pass
        else:
            new_paths.append(path)
            for i in range(2):
                data_path = path + "/freqChannel{}.txt".format(i+1)
                if os.path.exists(data_path):
                    mean, u = calc_uncer_from_file(data_path)
                    if u is None:
                        pass
                    else:
                        with open(data_path, "a")as file:
                            file.write("\n")
                            file.write("{} +- {}".format(mean, u))

def calc_uncer_from_file(path):
    with open(path, "r") as file:
        data = file.readlines()
    data_float = []
    for dp in data:
        data_float.append(float(dp.replace("\n", "")))
    sum = 0.0
    for dp in data_float:
        sum = sum + dp
    mean =sum / len(data_float)

    sigma = 0.0
    for dp in data_float:
        sigma = (dp - mean)**2
    if len(data_float) != 1 and len(data_float) != 0:
        sigma = sigma/(len(data_float)-1)/len(data_float)
        sigma = np.sqrt(sigma)
    else:
        sigma = None

    return mean, sigma


def FFT(data, step, Tmax=None):
    if Tmax is None: Tmax = 1e-7  # 1/fs
    if step < Tmax:
        maxlength = int(step / Tmax * len(data))
        length = len(data)
        data = downSample(data=data, maxlength=maxlength)
        step = int(length / maxlength) * step
    N = len(data)
    fft = np.fft.fft(data, N)
    data = None
    mag = abs(fft)*2/N #*2/N
    fft = None
    freq = np.fft.fftfreq(N, step)
    return mag, freq

def binarySearch(data, key):
    if len(data) == 1:
        return 0
    else:
        i = int(len(data)/2)-1
        if data[i] == key:
            return i
        elif data[i] > key:
            return binarySearch(data[0:i], key)
        elif data[i] < key:
            return i + binarySearch(data[i:-1], key)












if __name__ == '__main__':
    #TODO titledict only works if titles are keys are unique for every path
    
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

    #Autoplot(curr_path, titledict = titledict)
    measureFreqAuto(path=curr_path)
    save = False

    freqLim = (3e6, 5e6)
    path = "/home/ole/StorageServer/Work/BachelorThesis/RedPitaya16BitUndMessaufbau/Messergebnisse/11.01.2021/90degPulse"
    xlim = None # (0.0, 0.00003)
    
    #Plotter(paths=[path], title ="Test", save=save, xlim=xlim)
    #FFTAutoPlot(path=path, title="Output RP", xlim=xlim, save=save, freqLim=freqLim)




