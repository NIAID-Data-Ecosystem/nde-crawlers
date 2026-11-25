from biothings.hub.dataload.dumper import KubernetesJobDumper


class VDJ_Dumper(KubernetesJobDumper):
    SRC_NAME = "vdj"
    SCHEDULE = "0 20 * * 6"  # Every Saturday at 8:00 PM
    SRC_URLS = [
        "k8s://in-cluster?namespace=nde-hub-nonprod&podtemplate={SRC_NAME}-crawler-template&path=/data/{SRC_NAME}_crawled/data.ndjson"
        # Previous su07 setup (DockerContainerDumper):
        # f'docker://su07?image=nde-crawlers-{SRC_NAME}-crawler&tag=latest&path=/data/{SRC_NAME}_crawled/data.ndjson&dump_command="/home/biothings/run-api-crawler.sh"&container_name={SRC_NAME}_dumper&keep_container=false'
    ]
