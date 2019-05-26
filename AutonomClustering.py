from Bio import SeqIO
import os
from Bio.Align.Applications import ClustalwCommandline
from leven import levenshtein
import matplotlib.pyplot as plt
from math import sqrt
# 1.Дозапсиь в фаста!!!
clustalw_exe = r"C:\Program Files (x86)\ClustalW2\clustalw2.exe"


# def hamming()

def createcenter(k):                                       # Функция создает центры из исходного файла
    x = []
    for seq_record in SeqIO.parse(name_of_file, "fasta"):
        x.append(seq_record)
    seq_len = len(x)
    for i in range(k):                                     # Выбираем из исходного файла элементы
        f = open("old_center%s.txt" % i, "w")              # и записываем в txt
        for index in x[int((i * seq_len) / k)].seq:
            f.write(index)
        f.close()
    x.clear()
    for i in range(k):                                     # Создаем txt для новых центров
        f = open("new_center%s.txt" % i, "w")
        f.close()


def new_createcenter(k):
    seq_dict = SeqIO.index(name_of_file, "fasta")
    dict_len = len(seq_dict)
    ids_list = []
    iteration = []
    for i in range(dict_len):
        iteration.append(0)
    for i in range(k):
        iteration[int((i * dict_len) / k)] = 1
    intent = 0
    for rec in SeqIO.parse(name_of_file, "fasta"):
        if iteration[intent]:
            ids_list.append(rec.id)
        intent += 1
    intent = 0
    for id in ids_list:                                     # Выбираем из исходного файла элементы
        f = open("old_center%s.txt" % intent, "w")              # и записываем в txt
        for index in seq_dict[id].seq:
            f.write(index)
        f.close()
        intent += 1
    seq_dict.close()
    for i in range(k):                                     # Создаем txt для новых центров
        f = open("new_center%s.txt" % i, "w")
        f.close()


def compare(k):                                            # Функция сравнивает предыдущие и новые центры
    need_to_continue = 0                                   # и в случаем их совпадения выдает 0
    for i in range(k):
        oldf = open("old_centerK.txt".replace("K", str(i)), "r")
        old = oldf.read()
        oldf.close()
        newf = open("new_centerK.txt".replace("K", str(i)), "r")
        new = newf.read()
        newf.close()
        if old != new:
            need_to_continue = 1
        # old.clear()
        # new.clear()
    return need_to_continue


def cluster_number(seq, k):                                # Возвращает номер кластера к которому принадлежит seq
    oldf = open("new_center0.txt", "r")
    old_min = oldf.read()
    oldf.close()
    min_dist = levenshtein(str(seq), str(old_min))
    number = 0
    for n in range(k):
        oldf = open("new_centerK.txt".replace("K", str(n)), "r")
        old = oldf.read()
        oldf.close()
        dist = levenshtein(str(seq), str(old))
        if dist <= min_dist:
            min_dist = dist
            number = n
    print("Number:%s" % number)
    return number


def find_diam(k):                                         # Находит диаметр кластеров и записывает в txt
    write_diam = ""
    for i in range(k):
        max_diam = 0
        for seq_record_i in SeqIO.parse("clusterI.fasta".replace("I", str(i)), "fasta"):
            for seq_record_j in SeqIO.parse("clusterI.fasta".replace("I", str(i)), "fasta"):
                diam = levenshtein(str(seq_record_i.seq), str(seq_record_j.seq))
                if diam >= max_diam:
                    max_diam = diam
        write_diam += ("Maximum distance on cluster %s is %s" % (i, max_diam) + "\n")
    for i in range(k):
        f = open("diamter%s_info.txt" % k, "w")
        for index in write_diam:
            f.write(index)
        f.close()


def centering(name, n):                                   # Создает новые центры кластеров
    first_record = next(SeqIO.parse(name, "fasta"))
    center_general = ""
    sequance_lenght = len(first_record.seq)  # Узнаем длину выражения
    for i in range(sequance_lenght):
        a = t = g = c = space = 0
        for center in SeqIO.parse(name, "fasta"):
            if center.seq[i] == "A":
                a += 1
            elif center.seq[i] == "T":
                t += 1
            elif center.seq[i] == "G":
                g += 1
            elif center.seq[i] == "C":
                c += 1
            else:
                space += 1
        maximum = max(a, t, g, c, space)          # ?Можно запилить рандом
        if a == maximum:
            center_general += "A"
        elif t == maximum:
            center_general += "T"
        elif g == maximum:
            center_general += "G"
        elif c == maximum:
            center_general += "C"
        else:
            center_general += "-"
    f = open("old_centerN.txt".replace("N", str(n)), "w")
    for index in center_general:
        f.write(index)
    f.close()


def clustering(k):                                  # ?Сделать бы динамически ввод имени файла
    print(os.getcwd())
    for i in range(k):
        os.rename("old_centerI.txt".replace("I", str(i)), "dull.txt")
        os.rename("new_centerI.txt".replace("I", str(i)), "old_centerI.txt".replace("I", str(i)))
        os.rename("dull.txt", "new_centerI.txt".replace("I", str(i)))
    xui = []                                        # Создаем список для хранения Seq из файла
    for seq_record in SeqIO.parse(name_of_file, "fasta"):  # Заполняем массив элементами
        xui.append(seq_record)
    seq_list = []
    for i in range(k):
        seq_list.append([])
    for i in range(len(xui)):                   # Проверяем какому их двух центров выражение ближе и относим к кластеру
        seq_list[cluster_number(xui[i].seq, k)].append(xui[i])
    xui.clear()                                     # Удаляем элементы их списка
    for i in range(k):
        SeqIO.write(seq_list[i], "clusterI.fasta".replace("I", str(i)), "fasta")
        seq_list[i].clear()

# Можнно брать случайные 100 элементов выравнивать их и создавать центр
# И ослабить проверку остановки, пусть у нас будет критерием не совпадение, а отличие расстояния
def alig(k):
    for i in range(k):
        cluster_len = SeqIO.to_dict(SeqIO.parse("cluster%s.fasta" % i, "fasta"))
        if len(cluster_len) > 0:
            clustalw_cline = ClustalwCommandline(clustalw_exe, output="fasta",
                                                 infile="clusterI.fasta".replace("I", str(i)))  # Выравниваем 1 кластер
            assert os.path.isfile(clustalw_exe), "Clustal W executable missing"
            stdout, stderr = clustalw_cline()
        else:
            print("Empty cluster")


def rename(k):                                # Переопределяет центры кластеров
    for i in range(k):
        name = "clusterI.fasta".replace("I", str(i))
        centering(name, i)


def new_clustering(k):
    for i in range(k):
        os.rename("old_centerI.txt".replace("I", str(i)), "dull.txt")
        os.rename("new_centerI.txt".replace("I", str(i)), "old_centerI.txt".replace("I", str(i)))
        os.rename("dull.txt", "new_centerI.txt".replace("I", str(i)))
    #ids_list = (rec.id for rec in SeqIO.parse("ls_orchid.fasta", "fasta"))
    ids_list = []
    for rec in SeqIO.parse(name_of_file, "fasta"):
        ids_list.append(rec.id)
    new_id_list = []
    for i in range(k):
        new_id_list.append([])
    orchid_dict = SeqIO.index(name_of_file, "fasta")
    for id in ids_list:
        new_id_list[cluster_number(orchid_dict[id].seq, k)].append(id)
    for i in range(k):
        records = (orchid_dict[id] for id in new_id_list[i])
        SeqIO.write(records, "clusterI.fasta".replace("I", str(i)), "fasta")
    orchid_dict.close()


name_of_file = "sorted.fasta"
kmeans = 10  # Количество кластеро

def visualization_lenght():
    seq_len = []
    for seq in SeqIO.parse("Test.fasta", "fasta"):
        seq_len.append(len(seq.seq))
    plt.plot(seq_len)
    plt.show()

#for number_of_clusters in range(2, int(sqrt(kmeans)) + 2):
#    new_createcenter(number_of_clusters)
#    while compare(number_of_clusters):
#        new_clustering(number_of_clusters)
#        alig(number_of_clusters)
#        rename(number_of_clusters)
#    find_diam(number_of_clusters)

visualization_lenght()