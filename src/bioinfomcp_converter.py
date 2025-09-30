import os
import subprocess
import json
import asyncio
import requests
import ast # this is just for check code syntax whether or not it is correct
import re
import pymupdf
from dotenv import load_dotenv

load_dotenv()


'''
# 0. Check whether that particular library is already installed (otherwise cannot run the help function)
# 1. Run the --help function
# 2. Let GPT-4 analyze the result from the help manual, attain the
    a) What tools are there
        i) For each tools, what is the input, output
        ii) the CLI command, error handling
        iii) capture output
        iv) check output whether it is correct
'''


class BioinfoMCP():
    def __init__(self, model="openai"):
        file = open('system_prompt.txt', 'r')
        self.sys_prompt = file.read()
        self.api_model_name = os.getenv('MODEL_NAME')

        if model == "azure":
            from openai import AzureOpenAI
            # Load Environment Keys
            api_subscription_key = os.getenv('AZURE_OPENAI_KEY')
            api_version = os.getenv('AZURE_OPENAI_API_VERSION')
            api_endpoint = os.getenv('AZURE_OPENAI_ENDPOINT')
            # initialize the Azure Client
            self.client = AzureOpenAI(
                api_version=api_version,
                azure_endpoint=api_endpoint,
                api_key=api_subscription_key,
            )

        elif model == "openai":
            from openai import OpenAI
            # Load Environment Keys
            openai_api_key = os.getenv('OPENAI_API_KEY')
            # initialize the OpenAI Client
            self.client = OpenAI(
                api_key=openai_api_key
            )
        print(f"Succesfully created a {self.api_model_name} model")
        
        
    
    def is_tool_available(self, tool_name):
        """Check whether that tool is installed or not"""
        try:
            vers_result = subprocess.run([tool_name, '--version'],
                                         capture_output=True, timeout=20, text=True)
            print(f"{tool_name} is installed")
            return True
        except:
            print(f"{tool_name} is not installed!")
            return False
    
    def extract_help_document(self, tool_name, manual, run_help_command=False):
        """Extract help text from tool"""
        if not run_help_command: # document provided by user, no need to run the help command
            manual = pymupdf.open(f'{manual}')
            manual_content = ""
            for page in manual:
                text = page.get_text()
                manual_content += text 
            return manual_content
            
        if self.is_tool_available(tool_name):
            try:
                result = subprocess.run([tool_name, manual],
                                         capture_output=True, timeout=30, text=True)
                return result.stdout + result.stderr
            except:
                return None

    def generate_prompt(self, tool_name, help_docs):
        """Generate the Prompt that later will be send to the OpenAI client"""
        prompt = f"""
Convert the following bioinformatics tool into an MCP tool definition.

Tool Name: {tool_name}
Help Document:
{help_docs}

Parse the Input parameters correctly, follow the MCP best practices, and provide the complete python code with the @mcp.tool() decorator.
"""
        return prompt
    
    def parse_mcpcode(self, gpt_response):
        code_block = re.findall('```python\n(.*?)\n```', gpt_response, re.DOTALL)
        if not code_block:
            return (0, "There is no python code found in the gpt_response", None)
    
        code = code_block[0]
        try:
            ast.parse(code)
        except SyntaxError as e:
            return (0, f"Code is not working with the following SyntaxError: {e}", code)
    
        # Check the @mcp.tool() decorator
        if '@mcp.tool' not in code:
            return (0, "Code is missing the @mcp.tool() decorator", code)

        return (1, None, code)
    
    def refine_after_feedback(self, tool_name, code, error_message):
        prompt = f"""
            The initial code block for {tool_name}:
            ```python
            {code}
            ```
            contains the following Error:
            {error_message}

            Please fix the code and ensure that 
            1. Cover every internal functions of the tool
            2. Has proper error handling
            3. Validates input parameters
            4. Returns structured output
            5. Follows MCP best practices

            Provide only the corrected python code.

        """
        response = self.client.chat.completions.create(
            message=[
                {
                    "role" : "system",
                    "content": [{"type": "text", "text": self.sys_prompt}],
                },
                {
                    "role": "user",
                    "content": [{"type": "text", "text": prompt}],
                }
            ],
            model=self.api_model_name,
            temperature=0.1
        )
        response_content = response.choices[0].message.content
        # output_file = open(f'./mcp_result/raw_{tool_name}', 'w')
        # output_file.write(response_content)
        return self.parse_mcpcode(response_content)

    def autogenerate_mcp_tool(self, tool_name, manual, run_help_command):
        """Autogenerate the MCP Tool using the extracted help documents"""
        help_docs = self.extract_help_document(tool_name, manual, run_help_command)
        # web_docs = self.fetch_web_docs(tool_name)
        prompt = self.generate_prompt(tool_name, help_docs)
        # Call the OpenAI Client to do Chat Completion
        response = self.client.chat.completions.create(
            messages=[
                {
                "role" : "system",
                "content": [{"type": "text", "text": self.sys_prompt}],
            },
            {
                "role" : "user",
                "content": [{"type": "text", "text": prompt}],
            }
            ],
            temperature=0.1,
            model=self.api_model_name
        )
        response_content = response.choices[0].message.content
        #print(response_content)
        # output_file = open(f'./mcp_result/raw_{tool_name}', 'w')
        # output_file.write(response.choices[0].message.content)

        return self.parse_mcpcode(response_content)
