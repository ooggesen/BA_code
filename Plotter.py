#Created 12.03.2020 by Ole Oggesen

import numpy as np
import matplotlib.pyplot as plt
import os
import array
import xml.etree.ElementTree as ET


def plot(xarray, yarray, title, xlabel, ylabel,  save, labelarray = None, xlim = None, ylim = None, savePath = None):
    '''Function for plotting the measurments on the NMR MRT
    y can contain multiple dataarrays,
    save is either True or False
    xlim is a touple with 2 values
    Needs a plots folder in same folder to save the plots'''
    plt.figure()
    for x, y in zip(xarray, yarray):
        plt.plot(x, y, linewidth=1)
    #plt.title(title)
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
        if savePath is None:
            savetitle = "./plots/" + title.replace(" ", "_") + ".pdf"
        else:
            savePath = savePath + "/plots"
            if not os.path.exists(savePath):
                os.mkdir(savePath)
            savetitle = savePath + "/" + title.replace(" ", "_") + ".pdf"
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
            raise FileNotFoundError("File not found. Only .bin and .csv datatypes are supported.")
    val = {}

    # parse for x values
    if isinstance(args, type([])):
        for arg in args:
            if "BaseUnit:" in arg:
                tmp = arg.split(":")
                assert tmp[1] == "V"
            if "HardwareXStart" in arg:
                tmp = arg.split(":")
                val["HardwareXStart"] = float(tmp[1])
            elif "HardwareXStop" in arg:
                tmp = arg.split(":")
                val["HardwareXStop"] = float(tmp[1])

        #parse for y values
        if ";" in data[0]:
            tmp = []
            for i, point in enumerate(data):
                data_split = point.split(";")
                while len(tmp) < len(data_split):
                    tmp.append([])
                for j, d in enumerate(data_split):
                    tmp[j].append(float(d))
            data = tmp
        else:
            for i, point in enumerate(data):
                data[i] = float(point)
            data = [data]

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
    ylim = (np.min(yarray[0])*1.1, np.max(yarray[0])*1.1)

    plot(xarray, yarray, title, xlabel, "voltage in V", save, labelarray=labelarray, xlim=xlim, savePath=os.path.dirname(path), ylim=ylim)

def findData(path, search_key=None):
    """Finds the path of data in a given folder and its subfolders. If search_key is given, only paths which contain the search_key are returned.
    Else it returns all Wfm.bin and Wfm.csv files and excludes its endings."""
    paths = []
    with os.scandir(path) as it:
        for entry in it:
            if entry.is_dir():
                try:
                    paths.extend(findData(entry.path, search_key=search_key))
                except:
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
    """
    Muliplies the data in dataElem with a rectangularlar window.
    :param dataElem: array like containing information
    :param xlim: touple with start and end point
    :param args:
    :return:
    """
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
    """
    Downsampling.
    :param data: array like
    :param maxlength: maximum length of data array like structure
    :param x: Same as data for two dimensional information
    :return:
    """
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

def measureFreq(path, xlim=None, freqLim = None):
    """Makes a FFT of the received data, and writes the dominant frequency to an .txt file."""
    data, args = parseData(path)

    for i, d in enumerate(data):
        if xlim is not None:
            d, t, step = RectWindow(dataElem=d, xlim=xlim, args=args)
        else:
            step = (args["HardwareXStop"] - args["HardwareXStart"]) / len(d)
        T0 = step * len(d)
        T0min = 1/50    #Defines the sampling rate in the frequency domain: fs = 1/T=min TODO set bigger value if final run
        if T0min > T0:
            d.extend(np.zeros(int(T0min / step) - int(T0 / step)))      #Zero padding
        mag, freq = FFT(data=d, step=step)
        freq = freq[0:int(len(freq)/2)]
        if freqLim is None:
            freqLim = (3e6, 5e6)
        idx_start = binarySearch(freq, freqLim[0])
        idx_stop = binarySearch(freq, freqLim[1])

        m_max = 0.0
        for f, m in zip(freq[idx_start:idx_stop], mag[idx_start:idx_stop]):
            if m>m_max:
                f_max = f
                m_max = m

        file_path = os.path.dirname(path)
        with open(file_path + "/freqChannel{}.txt".format(i+1), "a") as file:
            file.write("{}\n".format(f_max))

def measureFreqAuto(path, xlim = None):
    """Measures the dominant frequency in the data stored in path address and saves it to a .txt file. Only frequencies between 3MhZ and 5 MhZ are evaluated.
    :param path: path and subfolders are searched for measuerment data
    :param xlim: (t_start, t_stop) Makes a rectangular window over given time period between t_start and t_stop
    :return: .txt file containing the frequency of highest amplitude between 3 MHz and 5 MHz"""
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
                    mean, u = calcUncerFromFile(data_path)
                    if u is None:
                        pass
                    else:
                        with open(data_path, "a")as file:
                            file.write("\n")
                            file.write("{} +- {}".format(mean, u))

def calcUncerFromFile(path):
    """
    Calculates the uncertainty and mean value.
    :param path: path to a .txt documents which contians the individual values
    :return: touple with mean and uncertainty information
    """
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


def FFT(data, step, Tmin=None):
    """
    Calculates a FFT. Capable of downsampling if step < Tmin .
    :param data: amplitude data
    :param step: time step
    :param Tmin: Minimum time step used for down sampling
    :return: touple containing magnitute and frequency information
    """
    if Tmin is None: Tmin = 1e-7  # 1/fs
    if step < Tmin:
        maxlength = int(step / Tmin * len(data))
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
    """
    Classical binary search.
    :param data: array like
    :param key: search key
    :return: index closest to key or key
    """
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

def measureAmp(path, Tmax = None, xlim = None):
    """Evaluates the amplitude of the measured data. Has a Threshold implemented, and only amplitudes over that threshold are detected.
    :param Tmax : declares the maximum time between two samples. If this value is bigger than in the data, the data is downsampled."""
    data, args = parseData(path)

    for i, y in enumerate(data):
        if xlim is not None:
            y, t, step =RectWindow(dataElem=y, xlim=xlim, args=args)
        else:
            step = (args["HardwareXStop"] - args["HardwareXStart"]) / len(y)
        if Tmax is None: Tmax = 1e-7  # 1/fs
        if step < Tmax:
            maxlength = int(step / Tmax * len(y))
            length = len(y)
            y = downSample(data=y, maxlength=maxlength)
            step = int(length / maxlength) * step
        dydt = np.gradient(y, step)

        amplitudes = []

        if i == 0:
            threshold = 0.04
        else:
            threshold = 2.5
        while len(amplitudes)<5:
            for j, dy in enumerate(dydt):
                if abs(dy) < 0.5 and abs(y[j]) > threshold:
                    amplitudes.append(abs(y[j]))
            threshold = threshold - 0.01



        amplitude = np.sum(amplitudes) / len(amplitudes)


        file_path = os.path.dirname(path)
        with open(file_path + "/ampChannel{}.txt".format(i+1), "a") as file:
            file.write("{}\n".format(amplitude))

def measureAmpAuto(path):
    """
    Measuers the amplitude of the measuerd data.
    :param path: given path and subfolders are searched for measured data
    :return: .txt file in same folder as measuerments containing the amplitude information and a mean and uncertainty value
    """
    paths = findData(path)
    for path in paths:  # Removes the stored data in .txt files if they exist
        path = os.path.dirname(path)
        for i in range(2):
            data = path + "/ampChannel{}.txt".format(i + 1)
            if os.path.exists(data):
                os.remove(data)
    for path in paths:
        measureAmp(path)
    new_paths = []
    for path in paths:
        path = os.path.dirname(path)
        if path in new_paths:
            pass
        else:
            new_paths.append(path)
            for i in range(2):
                data_path = path + "/ampChannel{}.txt".format(i + 1)
                if os.path.exists(data_path):
                    mean, u = calcUncerFromFile(data_path)
                    if u is None:
                        pass
                    else:
                        with open(data_path, "a")as file:
                            file.write("\n")
                            file.write("{} +- {}".format(mean, u))

def readValuesFromTxtAuto(path):
    """
    Reads all information from .txt files which were stored by the measuerFreqAuto or measureAmpAuto functions.
    :param path: given path and subfolders are searched for .txt files
    :return: A tree like structure containing all information from .txt file in location from path or subfolders
    """
    paths = findData(path=path, search_key=".txt")
    return_dict = {}
    for path in paths:
        data = readValuesFromTxt(path)
        path = path.split("/")
        path = path[1:]
        tmp = {path[-1] : data}
        for i in reversed(path[0:-1]):
            tmp = {i : tmp}
        return_dict = merge_dicts(tmp, return_dict)
    return return_dict

def merge_dicts(dict1, dict2):
    """
    Merges two dicts
    :param dict1:
    :param dict2:
    :return: The merged dictionary
    """
    if isinstance(dict1, type({})) and isinstance(dict2,  type({})):
        for key1 in dict1.keys():
            if key1 in dict2.keys():
                dict2[key1] = merge_dicts(dict1[key1], dict2[key1])
            else:
                dict2[key1] = dict1[key1]
        return dict2
    elif not isinstance(dict1, type({})) or not isinstance(dict2, type({})):
        return [dict1, dict2]








def readValuesFromTxt(path):
    """
    Reads data from .txt files which was evaluated by measureFreqAuto or measureAmpAuto functions.
    :param path: Path of the .txt file"
    :return: touple element containing the mean value and uncertainty value
    """
    if os.path.exists(path):
        with open(path, "r") as file:
            data = file.readlines()
            if len(data) == 0:
                raise AssertionError("There is no data contained in given path: {}".format(path))
        for i, point in enumerate(data):
            data[i] = point.replace("\n", "")
        data = data[-1]
        data = data.split("+-")
        for i, d in enumerate(data):
            data[i] = float(d)
        return data
    else:
        raise AssertionError("The file described by the path does not exist.")

def evaluate(path):
    """Extracts frequency and amplitude information of measurement data and returns it in a tree like structure.
    :param path = path to folder in which/ or in which subfolders, measuerments are stored"""
    measureFreqAuto(path=path)
    measureAmpAuto(path=path)
    return readValuesFromTxtAuto(path=path)

def plotMeasurementsWithUncertainty(data, title=""):
    """Reads information from the data and plots it. Plots are stored in "./plots" relative to the Plotter.py file.
    :param data = Output of readValuesFromTxt(), Tree like structure containing the data
    :return: plots in ./plots folder relative to Plotter.py"""
    if isinstance(data, type({})):
        flarmor = 4.3576
        ampMust = []
        teList = []
        rf90durationList =[]
        teFreqChannel1 = []
        durFreqChannel1 = []
        ampChannel1 = []
        for key in data.keys():
            if "mplitude" in key:
                if "Amplitude" in key:
                    amp = key.replace("Amplitude", "")
                else:
                    amp = key.replace("amplitude", "")
                amp = float(amp.replace("mV", ""))
                ampMust.append(amp)
                tmp = data[key]["ampChannel1.txt"]
                for i, t in enumerate(tmp):
                    tmp[i] = t*1000
                ampChannel1.append(tmp)
            elif "te" in key:
                te = key.replace("te", "")
                teList.append(float(te.replace("ms", "")))
                tmp = data[key]["freqChannel1.txt"]
                for i, t in enumerate(tmp):
                    if i == 0:
                        tmp[i] = t*1e-3 - flarmor*1e3
                    else:
                        tmp[i] = t*1e-3
                teFreqChannel1.append(tmp)
            elif "rf90duration" in key:
                rf90duration = key.replace("rf90duration", "")
                rf90durationList.append(float(rf90duration.replace("us", "")))
                tmp = data[key]["freqChannel1.txt"]
                for i, t in enumerate(tmp):
                    if i == 0:
                        tmp[i] = t*1e-3 - flarmor*1e3
                    else:
                        tmp[i] = t*1e-3
                durFreqChannel1.append(tmp)
            else:
                plotMeasurementsWithUncertainty(data=data[key], title=key)
        plotWithUncertainties(ampMust, ampChannel1, title + "Amplitude")
        plotWithUncertainties(teList, teFreqChannel1, title + "Te")
        plotWithUncertainties(rf90durationList, durFreqChannel1, title + "Rf90Duration")
    else:
        pass




def plotWithUncertainties(X, Y, title, save=True):
    """Plots the data in X and Y with errorbars. Each value of Y has two entires. The fist is the measured value. The second the uncertainty.
    X = data on x axis
    Y = data on y axis with uncertainty
    title = used for saving the figure
    save = if True saves the figure, if False shows the figure"""
    if X == [] or Y == []:
        pass
    else:
        for x, y in zip(X, Y):
            plt.errorbar(x, y[0], yerr=y[1], fmt=".k")
        if "Amplitude" in title:
            plt.xlabel("amplitude in mV")
            plt.ylabel("output in mV")
        elif "Te" in title:
            plt.xlabel("Te in ms")
            plt.ylabel("Delta f in kHz")
        elif "Rf90Duration" in title:
            plt.xlabel("rf90duration in us")
            plt.ylabel("Delta f in kHz")
        plt.grid()
        if save == True:
            savetitle = "./plots/" + title.replace(" ", "_") + ".pdf"
            plt.savefig(savetitle)
            plt.close("all")
        else:
            plt.show()



















if __name__ == '__main__':
    #TODO titledict only works if titles are keys are unique for every path
    
    titledict = {}

    curr_path = os.path.dirname(__file__)   #returns the current path of the python skript
    curr_path = os.path.dirname(curr_path)  # ".."
    #curr_path = curr_path + "/07042021RFAmplifier"
    Autoplot(curr_path, titledict = titledict)
    #measureFreqAuto(path=curr_path)
    #measureAmpAuto(path=curr_path)
    data = readValuesFromTxtAuto(path=curr_path)
    plotMeasurementsWithUncertainty(data=data)

    save = False

    freqLim = (3e6, 5e6)
    path = "/home/ole/StorageServer/Work/BachelorThesis/RedPitaya16BitUndMessaufbau/Messergebnisse/11.01.2021/90degPulse"
    xlim = None # (0.0, 0.00003)
    
    #Plotter(paths=[path], title ="Test", save=save, xlim=xlim)
    #FFTAutoPlot(path=path, title="Output RP", xlim=xlim, save=save, freqLim=freqLim)




