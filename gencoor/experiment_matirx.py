
class ExpMatrix:
    """Experimental Matrix is a class for loading, storing and \
       manipulating the files from different sources and format."""

    def __init__(self):
        self.names = []
        self.types = {}
        self.files = {}
        # self.labels = {}
        self.headers = []
        self.additional_columns = {}
        self.tags = {}

    def read(self, path):

        fixed_headers = ["name", "type", "file"]
        fixed_idx = {}
        with open(path) as f:
            for i, line in enumerate(f):
                if i == 0:
                    l = line.strip().split()
                    for h in fixed_headers:
                        if h in l:
                            fixed_idx[h] = l.index(h)
                    self.headers = [h for h in l if h not in fixed_headers]
                    if self.headers:
                        for h in self.headers:
                            self.additional_columns[h] = {}
                elif line.startswith("#"):
                    continue
                else:
                    l = line.strip().split()
                    if len(l) >= 4:
                        new_name = l[fixed_idx["name"]]
                        self.names.append(new_name)
                        self.types[new_name] = l[fixed_idx["type"]]
                        self.files[new_name] = l[fixed_idx["file"]]
                        # self.labels[new_name] = l[fixed_idx["label"]]
                        ll = list(set(range(len(l))) - set(fixed_idx.values()))
                        # print(range(len(l)))
                        # print()


                        if self.headers:
                            for j, h in enumerate(self.headers):
                                self.additional_columns[h][new_name] = l[ll[j]]
                    else:
                        continue

        self.extend_entries()
        self.sum_up_tags()

    def sum_up_tags(self):
        for n in self.names:
            self.tags[n] = []
            if self.headers:
                for h in self.headers:
                    self.tags[n].append(self.additional_columns[h][n])
            self.tags[n].append(self.types[n])
            # self.tags[n].append(self.labels[n])
            self.tags[n] = set(self.tags[n])

    def remove_a_name(self, name):
        self.names.remove(name)
        del self.types[name]
        del self.files[name]
        # del self.labels[name]
        if self.headers:
            for h in self.additional_columns.keys():
                del self.additional_columns[h][name]

    def duplicate_entry(self, name, new_name):
        self.names.append(new_name)
        self.types[new_name] = self.types[name]
        self.files[new_name] = self.files[name]
        # self.labels[new_name] = self.labels[name]
        if self.headers:
            for h in self.additional_columns.keys():
                self.additional_columns[h][new_name] = self.additional_columns[h][name]

    def extend_entries(self):
        def split_a_tag():
            for h in self.headers:
                for name, tag in self.additional_columns[h].items():
                    if "," in tag:
                        new_labels = tag.split(",")
                        combined_labels = [name+"_"+ s for s in new_labels]
                        for i, ll in enumerate(combined_labels):
                            self.duplicate_entry(name, ll)
                            self.additional_columns[h][ll] = new_labels[i]
                        self.remove_a_name(name)
                        return False
            return True

        def complete_all_tag():
            for h in self.headers:
                for name, tag in self.additional_columns[h].items():
                    if tag == ".":
                        tags = self.get_all_tags(h)
                        combined_labels = [name + "_" + s for s in tags]
                        for i, t in enumerate(combined_labels):
                            self.duplicate_entry(name, t)
                            self.additional_columns[h][t] = tags[i]
                        self.remove_a_name(name)
                        return False
            return True

        if self.headers:
            reach_end = False
            while not reach_end:
                reach_end = split_a_tag()
            reach_end = False
            while not reach_end:
                reach_end = complete_all_tag()

    def get_regions(self):
        res = []
        for name, type in self.types.items():
            if type == "regions":
                res.append(name)
        return res

    def get_signals(self):
        res = []
        for name, type in self.types.items():
            if type == "signals":
                res.append(name)
        return res

    def print(self):
        print("####################################")
        print("\t".join(["name", "type", "file"] + self.headers))
        for name in self.names:
            print("\t".join([name, self.types[name], self.files[name]] +
                            [self.additional_columns[h][name] for h in self.headers]))
        print("####################################")

    def get_all_tags(self, header):
        res = []
        if header in self.headers:
            for name in self.names:
                tag = self.additional_columns[header][name]
                if tag != "." and tag not in res:
                    res.append(tag)
        elif header == "regions":
            res = self.get_regions()
        elif header == "signals":
            res = self.get_signals()
        else:
            res = [""]
        return res

    def filter_by_tags(self, tags):
        res = []
        cue = [t for t in tags if t != ""]
        # print("cue: "+ " ".join(cue))
        # print(self.tags)
        for name in self.names:
            if set(cue) <= set(self.tags[name]):
                res.append(name)
        return res

    def get_file(self, name):
        return self.files[name]
    #
    # def get_label(self, name):
    #     return self.labels[name]




