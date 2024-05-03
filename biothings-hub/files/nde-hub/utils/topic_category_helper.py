import json


def add_topic_category(docs, source_name):
    """
    Adds 'topicCategory' to documents from a JSON file if the category list is not empty.
    Creates a file path from 'source_name', loading mappings of document IDs to topic lists.

    :param docs: List of dictionaries, each representing a document.
    :param source_name: String specifying the path to the JSON file containing topic categories.
    :return: List of documents, potentially updated with 'topicCategory'.
    """
    file_path = f"hub/dataload/sources/{source_name}/{source_name}.json"
    with open(file_path, 'r') as file:
        topic_categories = json.load(file)

    topic_dict = {list(item.keys())[0]: list(item.values())[0] for item in topic_categories}

    updated_docs = []
    for doc in docs:
        doc_id = doc['_id']
        topics = topic_dict.get(doc_id, [])
        if topics:
            doc['topicCategory'] = [
                {
                    "name": topic,
                    "curatedBy": {
                        "name": "gpt-3.5-turbo",
                        "url": "https://openai.com/index/chatgpt"
                    }
                } for topic in topics
            ]
        updated_docs.append(doc)
    return updated_docs
