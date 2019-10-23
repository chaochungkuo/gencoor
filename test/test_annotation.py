from gencoor.annotations import Annotation

#################################################################
#### Tests on Annotation
#################################################################

if __name__ == '__main__':
    ann = Annotation(name="test", genome="hg19")
    print(len(ann))